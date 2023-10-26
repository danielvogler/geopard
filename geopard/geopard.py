"""Geopard module."""

import logging
from csv import DictReader
from datetime import datetime
from math import asin, cos, radians, sin, sqrt
from typing import TypeVar

import numpy as np
from gpxpy import parse  # pip3 install gpxpy
from matplotlib import pyplot as plt
from scipy.interpolate import splev, splprep
from scipy.spatial.distance import cdist
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


class GeopardException(Exception):
    """Exception class."""


class GeopardResponse:
    """Response class."""

    time = None
    dtw = 0
    start_point = None
    end_point = None
    match_flag = 0
    error = 0

    def __init__(self, time, dtw, start_point, end_point, match_flag, error=None):
        """init."""
        self.time = time
        self.dtw = dtw
        self.start_point = start_point
        self.end_point = end_point
        self.match_flag = match_flag
        self.error = error

    def time(self):
        """Time."""
        return self.time

    def dtw(self):
        """Dynamic time warping metric."""
        return self.dtw

    def start_point(self):
        """Start point."""
        return self.start_point

    def end_point(self):
        """End point."""
        return self.end_point

    def match_flag(self):
        """Match flag."""
        return self.match_flag

    def error(self):
        """Error."""
        return self.error

    def is_success(self):
        """Success flag."""
        return self.match_flag > 0


class Geopard:
    """Main class."""

    def dtw_match(
        self,
        gold_name: str,
        activity_name: str,
        min_trkps: int = 50,
        radius: float = 7.0,
        dtw_threshold: float = 0.2,
        dtw_margin_range: float = 1.5,
        start_region=None,
        finish_region=None,
    ) -> TypeVar:
        """Prepare gpx tracks for dtw matching.

        Find all suitable start/end point combinations. determine shortest segment match.

        Args:
            gold_name (str): Name of gold segment
            activity_name (str): Name of activity segment
            min_trkps (int, optional): Minimum number of GPS points. Defaults to 50.
            radius (float, optional): Radius in meters. Defaults to 7.
            dtw_threshold (float, optional): DTW threshold. Defaults to 0.2.
            dtw_margin_range (float, optional): DTW margin range. Defaults to 1.5.
            start_region (bool, optional): _description_. Defaults to None.
            finish_region (_type_, optional): _description_. Defaults to None.

        Raises:
            GeopardException: Potential exception

        Returns:
            GeopardResponse: succesful or unsuccesful response
        """
        logging.info("Starting DTW match")
        start_time = datetime.now()

        ### load gold standard/baseline segment
        gold = self.gpx_loading(gold_name)
        ### interpolate gold data
        gold_interpolated = self.interpolate(gold)

        ### load activity data to be edited
        trkps = self.gpx_loading(activity_name)

        ### check if gpx track [lat, lon, ele, time] is complete. Missing elevation data is accepted
        for k in [0, 1, 3]:
            none_idx = [i for i in range(len(trkps[k])) if trkps[k, i] == None]

            ### error if required gpx track information is missing
            if none_idx:
                raise GeopardException(
                    "GPX file corrupted. Check GPX points (E.g. {})".format(sorted(set(none_idx))[:3])
                )

        ### search activity trackpoints matching with gold start/end
        try:
            ### crop activity data to segment length
            gpx_cropped = self.gpx_track_crop(gold, trkps, radius, start_region, finish_region)

            ### find potential start/end trackpoints
            nn_start, nn_start_idx = self.nearest_neighbours(gpx_cropped, gold[:4, 0], radius, start_region)
            nn_finish, nn_finish_idx = self.nearest_neighbours(gpx_cropped, gold[:4, -1], radius, finish_region)

        except GeopardException as e:
            return GeopardResponse(None, None, None, None, -2, str(e))

        ### tested combinations of start/end points
        combinations_to_test = 0

        ### initialize segment combinations
        seg_time = []
        seg_start_idx = []
        seg_finish_idx = []

        ### possible start point combinations
        for i in nn_start_idx:
            ### possible end point combinations
            for j in nn_finish_idx:
                ### start needs to happen before finish and include min_trkps in between
                if i < j - min_trkps:
                    combinations_to_test += 1

                    ### segment_time
                    seg_time.append(gpx_cropped[3, j] - gpx_cropped[3, i])
                    seg_start_idx.append(i)
                    seg_finish_idx.append(j)

        ### compile information of individual times and array indices
        seg_info = np.asarray([seg_time, seg_start_idx, seg_finish_idx])

        ### find indices of shortest runs to loop through in ascending order
        seg_sort = [seg_info[:, idx] for idx, value in sorted(enumerate(seg_info[0]), key=lambda x: x[1])]

        ### initialize segment check parameters
        s = -1
        final_dtw = 1e6
        match_flag = -1
        dtw_threshold_soft = dtw_threshold * dtw_margin_range
        dtw = None

        ### compute dtw in ascending order
        ### return if shortest activity satisfies dtw requirement
        ### exit if all combinations are tested
        while final_dtw > dtw_threshold and combinations_to_test > s + 1:
            ### counter
            s += 1

            ### compute DTW between gold and activity
            dtw, delta_time = self.dtw_computation(
                gpx_cropped[:, seg_sort[s][1] : seg_sort[s][2] + 1], gold_interpolated
            )

            logging.warning("DTW (y): %2.5f / T [s]: %s (%s/%s)" % (dtw, delta_time, s, combinations_to_test))

            ### update dtw information if threshold is crossed
            if dtw <= dtw_threshold:
                ### update final time and dtw
                final_time = delta_time
                final_dtw = dtw
                final_start_point = gpx_cropped[:, seg_sort[s][1]]
                final_end_point = gpx_cropped[:, seg_sort[s][2]]
                match_flag = 2

            ### check if soft dtw threshold is observed.
            ### only save dtw value and time for shortest (first) match in grey zone
            ### shortest grey zone match is overwritten if hard dtw threshold is crossed
            elif dtw_threshold < dtw <= dtw_threshold_soft and match_flag < 0:
                ### update final time and dtw
                final_time = delta_time
                final_dtw = dtw
                final_start_point = gpx_cropped[:, seg_sort[s][1]]
                final_end_point = gpx_cropped[:, seg_sort[s][2]]
                match_flag = 1

        logging.info("----- Finished DTW segment match -----")

        logging.info("Total combinations to test: %s" % (combinations_to_test))
        logging.info("Total combinations tested: %s" % (s + 1))

        logging.info("Total execution time: %s" % str(datetime.now() - start_time))

        if s > -1:
            logging.info("Execution time per combination: %s" % str((datetime.now() - start_time) / (s + 1)))

        logging.info("Match flag [-]: %s" % (match_flag))

        logging.info("Checking if segment match was achieved.")
        if match_flag > 0:
            logging.info("Final T [s]: %s" % (final_time))
            logging.info("Final DTW (y): %2.5f" % (final_dtw))

            return GeopardResponse(final_time, final_dtw, final_start_point, final_end_point, match_flag)

        else:
            logging.warning("No segment match found")

            return GeopardResponse(None, dtw, None, None, match_flag)

    def gpx_loading(self, file_name: str) -> TypeVar:
        """Load gpx file and extract information."""
        logging.info(f"Load gpx file: {file_name}")

        ### load file
        gpx_data_open = open(file_name)
        gpx_data = parse(gpx_data_open)

        ### initialize lat/lon
        trkp_lat = []
        trkp_lon = []
        trkp_ele = []
        trkp_time = []

        logging.info(f"\tLoading relevant gpx data")
        for track in gpx_data.tracks:
            for segment in track.segments:
                for i in range(0, len(segment.points)):
                    trkp_point = segment.points[i]
                    trkp_lat.append(trkp_point.latitude)
                    trkp_lon.append(trkp_point.longitude)
                    trkp_ele.append(trkp_point.elevation)
                    trkp_time.append(trkp_point.time)

        logging.info(f"\tChecking for duplicate lat/lon trackpoints")
        duplicates_lat = [
            i for i in range(1, len(trkp_lat) - 1) if trkp_lat[i] == trkp_lat[i - 1] and trkp_lat[i] == trkp_lat[i + 1]
        ]
        duplicates_lon = [
            i for i in range(1, len(trkp_lon) - 1) if trkp_lon[i] == trkp_lon[i - 1] and trkp_lon[i] == trkp_lon[i + 1]
        ]

        ### duplicate values for both lat and lon
        duplicates_lat_lon = sorted(list(set(duplicates_lat) & set(duplicates_lon)))

        logging.info(f"\tRemoving redundant trackpoints")
        trkp_lat = np.delete(trkp_lat, duplicates_lat_lon, 0)
        trkp_lon = np.delete(trkp_lon, duplicates_lat_lon, 0)
        trkp_ele = np.delete(trkp_ele, duplicates_lat_lon, 0)
        trkp_time = np.delete(trkp_time, duplicates_lat_lon, 0)

        ### merge lat, lon, ele, time into one array
        trkps = np.asarray([trkp_lat, trkp_lon, trkp_ele, trkp_time])

        return trkps

    ### filter gpx data occurring between two points
    def gpx_track_crop(
        self,
        gold: TypeVar,
        gpx_data: TypeVar,
        radius: float = None,
        start_region: TypeVar = None,
        finish_region: TypeVar = None,
    ) -> TypeVar:
        """Cropping GPX track.

        Args:
            gold (TypeVar): gold segment
            gpx_data (TypeVar): gpx track data
            radius (float, optional): Radius. Defaults to None.
            start_region (TypeVar, optional): Start region of track. Defaults to None.
            finish_region (TypeVar, optional): Finish region of track. Defaults to None.

        Returns:
            TypeVar: _description_
        """
        logging.info("Cropping GPX track.")

        ### find possible start/end trackpoints
        nn_start, nn_start_idx = self.nearest_neighbours(gpx_data, gold[:4, 0], radius, start_region)
        nn_finish, nn_finish_idx = self.nearest_neighbours(gpx_data, gold[:4, -1], radius, finish_region)

        ### determine time range for gpx track
        nn_start_earliest = min(nn_start[3])
        nn_finish_latest = max(nn_finish[3])

        ### delete trackpoints outside relevant time range
        indices = [i for i in range(len(gpx_data[3])) if nn_start_earliest <= gpx_data[3, i] <= nn_finish_latest]
        cr_lat = [gpx_data[0, i] for i in indices]
        cr_lon = [gpx_data[1, i] for i in indices]
        cr_ele = [gpx_data[2, i] for i in indices]
        cr_time = [gpx_data[3, i] for i in indices]
        gpx_cropped = np.asarray([cr_lat, cr_lon, cr_ele, cr_time])

        return gpx_cropped

    ### find nearest neighbouring points of input point
    def nearest_neighbours(self, gpx_data: TypeVar, centroid=list, radius: float = None, region: TypeVar = None):
        """Find nearest neighbors of GPS point.

        Args:
            gpx_data (TypeVar): _description_
            centroid (TypeVar, optional): Centroid around which neighbors should be found. Defaults to None.
            radius (float, optional): Radius. Defaults to None.
            region (TypeVar, optional): Region object. Defaults to None.

        Raises:
            GeopardException: _description_

        Returns:
            list: list with nearest neighbors
        """
        logging.info("Finding nearest neighbours ...")

        ### polygon region is given to find NN
        if region is not None:
            ### Convert track coordinates to points
            points = [Point(reversed(i)) for i in gpx_data[:4, :][:2].T]

            idx = []
            ### go through all gpx track points
            for i in range(len(points)):
                ### check if track points are within region polygon
                if region.contains(points[i]):
                    idx.append(int(i))

            logging.info("{} NN within region polygon: {}".format(len(idx), region))

        ### no polygon region given - default to circle around start
        else:
            ### distance (m) of all points to centroid
            distance = [self.spheroid_point_distance(i, centroid[:2]) for i in gpx_data[:4, :][:2].T]

            ### points within radius distance of centroid
            idx = [int(i) for i, x in enumerate(distance) if x < radius]

            logging.info("{} NN within radius of {}m near centroid: {}".format(len(idx), radius, centroid[:2]))

        ### filter out consecutive candidates
        nums = sorted(set(idx))
        gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s + 1 < e]
        edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
        idx = list(edges)

        if not idx:
            raise GeopardException("No trackpoints found near centroid")

        ### lat, lon, ele, time, distance of all nearest neighbours
        lat = [gpx_data[0, i] for i in idx]
        lon = [gpx_data[1, i] for i in idx]
        ele = [gpx_data[2, i] for i in idx]
        time = [gpx_data[3, i] for i in idx]
        # dis   = [distance[i] for i in idx]
        nn = np.asarray([lat, lon, ele, time])

        return nn, idx

    ### determine great-circle distance of two coordinate points on a sphere
    def spheroid_point_distance(self, coordinates_1: list, coordinates_2: list) -> float:
        """Compute distance of two points on spheroid.

        Args:
            coordinates_1 (list): Point 1
            coordinates_2 (list): Point 2

        Returns:
            float: Distance between both points
        """
        logging.debug(f"\t Compute distance between two points on spheroid.")

        avg_earth_radius = 6371008.7714  # [m] https://doi.org/10.1007/s001900050278

        lat_1 = radians(coordinates_1[0])
        lat_2 = radians(coordinates_2[0])
        lon_1 = radians(coordinates_1[1])
        lon_2 = radians(coordinates_2[1])

        delta_lat = (lat_2 - lat_1) / 2
        delta_lon = (lon_2 - lon_1) / 2

        ### from haversine formlua: https://en.wikipedia.org/wiki/Haversine_formula
        distance = (
            2 * avg_earth_radius * asin(sqrt(sin(delta_lat) ** 2 + cos(lat_1) * cos(lat_2) * sin(delta_lon) ** 2))
        )

        return distance

    ### accumulated cost matrix
    def acm(self, reference, query, distance_metric: str = "euclidean") -> TypeVar:
        """Accumulated cost matrix calculation.

        According to:
        (1) Müller, Meinard. Information retrieval for music and motion. Vol. 2.
            Heidelberg: Springer, 2007. https://doi.org/10.1007/978-3-540-74048-3
        """
        logging.debug(f"\t Compute ccumulated cost matrix")

        # cm = self.cost_matrix(reference, query)
        cm = cdist(reference, query, metric=distance_metric)

        ### sequence lengths
        N, M = cm.shape

        ### initialize
        acm = np.zeros([N, M])

        ### boundary condition 1
        acm[0, 0] = cm[0, 0]

        ### From (1) Theorem 4.3
        ### D(n, 1) = \sum_{k=1}^n c(x_k , y_1 ) for n ∈ [1 : N ],
        for n in range(1, N):
            acm[n, 0] = acm[n - 1, 0] + cm[n, 0]

        ### D(1, m) = \sum_{k=1}^n c(x_1 , y_k ) for m ∈ [1 : M ] and
        for m in range(1, M):
            acm[0, m] = acm[0, m - 1] + cm[0, m]

        ### for 1 < n ≤ N and 1 < m ≤ M .
        ### D(n, m) = min{D(n − 1, m − 1), D(n − 1, m), D(n, m − 1)} + c(x_n , y_m )
        for n in range(1, N):
            for m in range(1, M):
                acm[n, m] = cm[n, m] + min(acm[n - 1, m], acm[n, m - 1], acm[n - 1, m - 1])

        return acm

    def dtw_computation(self, gpx_data: TypeVar, gold: TypeVar) -> list:
        """Compute dynamic time warping between two gpx-track curves.

        Args:
            gpx_data (TypeVar): Activity track
            gold (TypeVar): Gold standard track

        Returns:
            list: DTW and time
        """
        logging.debug("Compute dynamic time warping between two gpx-track curves.")

        ### time difference of start to finish
        delta_time = gpx_data[3][-1] - gpx_data[3][0]

        ### interpolate activity
        gpx_data_interpolated = self.interpolate(gpx_data)

        ### compute dynamic time warping
        dtw = self.acm(gpx_data_interpolated, gold)[-1, -1]

        return dtw, delta_time

    def interpolate(self, gpx_data: TypeVar) -> TypeVar:
        """Interpolate gpx data along track.

        Args:
            gpx_data (TypeVar): GPX track

        Returns:
            TypeVar: interpolated track
        """
        logging.debug("Interpolate gpx data along track.")

        interpolation_points = 1000

        ### add random noise to avoid duplicates
        data = np.array(
            [
                (x + np.random.uniform(low=-1e-7, high=1e-7), y + np.random.uniform(low=-1e-7, high=1e-7))
                for x, y in zip(gpx_data[0][:], gpx_data[1][:])
            ]
        )

        ### function that interpolates original data
        tck, u = splprep(data.T, u=None, s=0.0, t=10, per=0)

        ###
        u_new = np.linspace(u.min(), u.max(), interpolation_points)

        ### fit interpolation function to equi-distant data
        lat, lon = splev(u_new, tck, der=0)

        ### only return necessary gpx data
        points = np.zeros((interpolation_points, 2))

        # points = np.asarray([lat, lon])
        points[:, 0] = lat
        points[:, 1] = lon

        return points

    def create_polygon(self, file_name: str) -> TypeVar:
        """Create polygon geometry from csv files.

        Args:
            file_name (str): file name containing polygon information

        Returns:
            TypeVar: polygon geometry
        """
        logging.debug("Create polygon")

        region = []
        ### read in lat/lon and conver to points
        with open(file_name, newline="") as f:
            reader = DictReader(f)
            for row in reader:
                region.append(Point([float(row["Longitude"]), float(row["Latitude"])]))

        ### convert points to polygon
        poly = Polygon([[p.x, p.y] for p in region])

        return poly

    ### plot gpx track
    def gpx_plot(self, fig: TypeVar, gpx_data: TypeVar, plot_info: str, marker_size: int = 500) -> TypeVar:
        """Plot GPX data.

        Args:
            fig (TypeVar): Figure object
            gpx_data (TypeVar): GPX data
            plot_info (str): Plotting info
            marker_size (int, optional): Marker size. Defaults to 500.

        Returns:
            TypeVar: Figure object
        """
        logging.debug("Plot gpx data")

        font_size = 30
        plt.rcParams.update({"font.size": font_size})
        plt.scatter(
            gpx_data[1, :], gpx_data[0, :], s=marker_size, marker=plot_info[1], c=plot_info[2], label=plot_info[0]
        )
        plt.legend(loc="upper left")
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")

        return fig
