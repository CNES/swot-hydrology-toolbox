"""
 This file is part of the SWOT Hydrology Toolbox
 Copyright (c) 2026 CNES
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY,
 read LICENSE.txt for further details.
"""

import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

print_debug = False


def plot_nb_obj_as_a_function_of_water_body_size(water_bodies_size_tab, water_body='river width',
                                                 begin_slices_list=None, end_slices_list=None):
    """
    Plot the number of objects (reaches or lakes) by slices of water body size (river width or lake size)
    :param water_bodies_size_tab: list of water body size (river width or lake size)
    :type water_bodies_size_tab: list
    :param water_body: water body type, among 'river width' or 'lake size' (default)
    :type water_body: str
    :param begin_slices_list: list of begin of slices of water body sizes
    :type begin_slices_list: 1D-array of float
    :param end_slices_list: list of end of slices of water body sizes
            if len(end_slices_list) == len(begin_slices_list) - 1, the last slice corresponds to > begin_slices_list[-1]
    :type end_slices_list: 1D-array of float
    :return:
    """
    # 1.1 - Prepare data and X-axis title
    nb_slices = len(begin_slices_list)
    _slices_title_list = []
    tab_nb_items = np.zeros(nb_slices)
    for s, cur_slice in enumerate(begin_slices_list):

        # Build list of indices corresponding to slices
        index_sup = (water_bodies_size_tab > int(cur_slice))
        if (s == nb_slices - 1) and (end_slices_list[s] == "+"):
            water_bodies_size_selection = index_sup
        else:
            # Select water_bodies_size_tab index where size value is <= end slice
            index_inf = (water_bodies_size_tab <= int(end_slices_list[s]))
            # Select water_bodies_size_tab index where both conditions are gathered
            water_bodies_size_selection = np.logical_and(index_sup, index_inf)
        tab_nb_items[s] = len(np.where(water_bodies_size_selection)[0])

        # Build text of water body size (for plotting X labels)
        if (s == nb_slices - 1) and (end_slices_list[s] == "+"):
            slice_title = '[>{0}{1}]'.format(int(np.sqrt(cur_slice)), '\u00B2')  # for example: >1000\u00B2
        else:
            slice_title = '[{0}{1}-{2:d}{3}]'.format(int(np.sqrt(cur_slice)), '\u00B2',
                                                     int(np.sqrt(end_slices_list[s])),
                                                     '\u00B2')  # for example: [50\u00B2-100\u00B2]
        # Add slice title to the list
        _slices_title_list.append(slice_title)
    nb_total = np.sum(tab_nb_items)
    if nb_total == 0:
        return None

    # 1.2 - Compute related array of x indexes
    ind_x = range(nb_slices)

    # 2 - Init figure
    fig = plt.figure(figsize=(10, 4))
    plt.ylim(0, np.max(tab_nb_items) * 1.1)

    # 3 - Plot histogram of nb of items vs water body size slice
    barlist = plt.bar(ind_x, tab_nb_items)
    plt.xticks(ind_x, _slices_title_list)
    # Define a different color for each water bodies size
    if len(_slices_title_list) == 4:
        colorbar = {0: 'darkviolet', 1: 'violet', 2: 'darkorange', 3: 'gold'}
    else:
        colorbar = {0: 'darkblue', 1: 'darkviolet', 2: 'violet', 3: 'darkorange', 4: 'gold', 5: 'darkkhaki'}
    for ind in ind_x:
        barlist[ind].set_color(colorbar[ind])

    # 4 - Add title and labels
    # 4.1 - Title
    if water_body == "river width":
        item_name = "reaches"
    else:
        item_name = "lakes"
    title = 'Number of %s as a function of %s (%d in total)' % (item_name, water_body, nb_total)
    plt.title(title)
    # 4.2 - Label X
    if water_body == 'river width':
        water_body_unit = 'm'
    else:
        water_body_unit = 'm' + '\u00B2'
    xlabel = '{0} ({1})'.format(water_body.capitalize(), water_body_unit)
    plt.xlabel(xlabel)
    # 4.3 - Label Y
    plt.ylabel('Number of {}'.format(item_name))
    # 4.4 - Nb and associated percentage
    for ind, value in enumerate(tab_nb_items):
        print_str = "{0:.0f}%".format(round(value / nb_total * 100.0))
        if max(tab_nb_items) <= 100:
            coeff = 0.1
        else:
            coeff = 1
        plt.text(ind - 0.1, value + coeff, print_str, color=colorbar[ind], fontweight='bold')

    return fig


class MeasurementError(object):
    """ Class for displaying measurement errors, computing one sigma error(s)"""

    def __init__(self, errors, metrix_name='area_total', water_bodies_size=None, xrange=None):
        """
        Initialize measurement errors

        Args:
            errors (list): list of errors
            metrix_name: metrix_name, among 'wse' or 'area_total' or 'area_geom'
            water_bodies_size (list): list of water body size (river witdh or lake size)
            xrange (list): crosstrack distance of centroid of water feature; one-to-one linked to errors
        """

        # Save self variables wrt input parameters
        self._errors = errors
        self._metrix_name = metrix_name
        self._water_bodies_size = water_bodies_size
        self._xrange = xrange

        # Initialization
        self._errors_one_sigma = []
        self._index_sliced = dict()
        self._slices_title_list = []
        self._req_th_legend = None
        # Depending on metric
        if self._metrix_name == 'wse':
            self.metrix_long_name = "water surface elevation"
            self.error_type = ''  # Absolute but don't write it in graphs
            self.err_req = 0.10  # Requirement level = 10cm for lakes surface > 1km²
            self.err_goal = 0.25  # Requirement level = 25cm  for lakes surface between 250²m² and 1km²
            self.unit = 'm'
            self.factor = 1
        else:
            self.metrix_long_name = "water surface area"
            self.error_type = 'Relative'
            self.err_req = 0.15 * 100  # Requirement level for river width > 100m or lake surface > 250²m²
            self.err_goal = 0.25 * 100  # Goal level  for river width between [50-100]m or lake surface [100²-250²]m²
            self.unit = '%'
            self.factor = 100

    def print_metrix_value(self, in_value):
        """
        Convert in_value into a string approximated number and its corresponding unit
        
        :param in_value: value to convert
        :type in_value: float
        
        :return: out_str_value = converted value
        :rtype: out_str_value = str
        """

        # 1 - Set the correct decimal approximation
        if self._metrix_name == 'wse':
            out_str_value = "%.2f" % float(in_value)
        else:
            out_str_value = "%.0f" % float(in_value)

        # 2 - Add the unit
        out_str_value += self.unit

        return out_str_value

    def compute_error_one_sigma(self, threshold=None, begin_slices_list=None,
                                end_slices_list=None, water_body='river width', which_perfos='req'):
        """
        1-sigma error: computes the 68.27th percentile of the empirical
        error distribution.

        :param threshold: defined threshold option to ignore abnormal error values
        :type threshold: float
        :param begin_slices_list: list of begin of slices of water body sizes
        :type begin_slices_list: 1D-array of float
        :param end_slices_list: list of end of slices of water body sizes
                if len(end_slices_list) == len(begin_slices_list) - 1, the last slice corresponds to > begin_slices_list[-1]
        :type end_slices_list: 1D-array of float
        :param water_body: water body type, among 'river width' and 'lake size'
        :type water_body: str
        :param which_perfos: 'req' or 'goal'
        :type which_perfos: str
        Set:
        self._errors_one_sigma: 1-sigma error
        1D-array of float, of size len(end_slices_list)+1; 1st value corresponds to the global 1sigma error
        self._slices_title_list : titles corresponding to list of slices
        """
        logger = logging.getLogger(self.__class__.__name__)

        abs_errors = abs(self._errors)

        # 1 - Init variables
        # Convert water_bodies_size list into array
        water_bodies_size_tab = np.asarray(self._water_bodies_size)
        # Selected index to compute globals 1sigma error
        ind_for_global_1 = set()  # indexes for rivers width < 100m and lake surface < (250m)² if metric_name='area'
        # indexes for lake surface between (250m)² and 1km² if metric_name='wse'
        ind_for_global_2 = set()  # indexes for rivers width > 100m and lake surface > (250m)² if metric_name='area'
        # indexes for lake surface > 1km² if metric_name='wse'
        ind_for_global_3 = set()  # indexes for lake surface > 1km² if metric_name='area'
        ind_for_global = set()  # set of indexes to compute global mean/median values
        # Nb slices
        nb_slices = len(begin_slices_list)
        for s, cur_slice in enumerate(begin_slices_list):  # Loop over slices
            # 2 - Consider criteria

            # 2.1 - Select indices corresponding to the current slice of water body size
            # Select water_bodies_size_tab index where size value is > begin slice
            index_sup = (water_bodies_size_tab > np.int32(cur_slice))
            if (s == nb_slices - 1) and (end_slices_list[s] == "+"):
                water_bodies_size_selection = index_sup
            else:
                # Select water_bodies_size_tab index where size value is <= end slice
                index_inf = (water_bodies_size_tab <= np.int32(end_slices_list[s]))
                # Select water_bodies_size_tab index where both conditions are gathered
                water_bodies_size_selection = np.logical_and(index_sup, index_inf)

            # 2.2 - Add criteria on the threshold (or not) and list indices of interest
            if threshold is None:
                # All errors values are considered
                index_selected = np.where(water_bodies_size_selection)
                if s == 0:
                    logger.info("NO threshold")
            else:
                # Furthermore, delete abnormal errors values using a threshold
                index_selected = np.where(np.logical_and(water_bodies_size_selection, (abs_errors < threshold)))
                if s == 0:
                    logger.info("threshold=%f" % threshold)
            # Save into self variable
            self._index_sliced[s] = index_selected
            # 2.3 - Save indices for global 1 sigma errors:
            #  |--- if which_perfos='goal' : exclude 1st slice [50-100]/[100²-250²] to compute global 1 sigma error
            #  |--- if which_perfos='req' : don't exclude 1st slice
            #  |--- if metrix_name='wse' (only for lakes) :
            #                               compute ind_for_global_1 for lake surfaces between 250²m² and 1km2
            #                               compute ind_for_global_2 for lake surfaces between > 1km2
            logger.info("Perfos is computed in relation to %s" % which_perfos)

            if self._metrix_name.startswith("area_"):
                if which_perfos == 'goal':
                    if s == 0:
                        # indexes for rivers width < 100m and lake surface < (250m)²
                        pinfo = "%s/%s - indexes for rivers width < 100m and lake surface < (250m)²" % (
                        self._metrix_name, which_perfos)
                        ind_for_global_1 = ind_for_global_1.union(list(index_selected[0]))
                    else:
                        # indexes for rivers width > 100m and lake surface > (250m)²
                        pinfo = "%s/%s - indexes for rivers width > 100m and lake surface > (250m)²" % (
                        self._metrix_name, which_perfos)
                        ind_for_global_2 = ind_for_global_2.union(list(index_selected[0]))
                    if s >= 3:
                        pinfo = "%s/%s - indexes for lake surface > 1km²" % (self._metrix_name, which_perfos)
                        ind_for_global_3 = ind_for_global_3.union(list(index_selected[0]))
                else:
                    # indexes for rivers width > 100m and lake surface > (250m)²
                    pinfo = "%s/%s - indexes for rivers width > 100m and lake surface > (250m)²" % (
                    self._metrix_name, which_perfos)
                    ind_for_global_2 = ind_for_global_2.union(list(index_selected[0]))
                    ind_for_global_1 = ind_for_global_2  # not use

                # set requirement threshold legend
                if water_body == 'lake size':
                    self._req_th_legend = '(250m){0}'.format('\u00B2')  # (250m)²
                else:  # for river width
                    self._req_th_legend = '100m'
            else:
                # self._metrix_name == 'wse'
                if which_perfos == 'goal':
                    s_lim_max = 3  # 1km²
                else:
                    s_lim_max = 2  # 1km²
                # ignore indexes < (250m)²
                if s >= 0 and s < s_lim_max:
                    # goal: indexes for lake surface between (100m)² and 1km²
                    # req: indexes for lake surface between (250m)² and 1km²
                    if which_perfos == 'goal':
                        pinfo = "%s/%s - indexes for lake surface between (100m)² and 1km²" % (
                        self._metrix_name, which_perfos)
                    else:
                        pinfo = "%s/%s - indexes for lake surface between (250m)² and 1km²" % (
                            self._metrix_name, which_perfos)
                    ind_for_global_1 = ind_for_global_1.union(list(index_selected[0]))
                else:
                    # indexes for lake surface > 1km²
                    pinfo = "%s/%s - indexes for lake surface > 1km²" % (self._metrix_name, which_perfos)
                    ind_for_global_2 = ind_for_global_2.union(list(index_selected[0]))
                self._req_th_legend = '1km{0}'.format('\u00B2')  # 1km²

            # 3 - Compute measurement error one sigma for slice s
            if len(abs_errors[index_selected]) == 0:
                error_one_sigma = 0  # No data found into this slice
            else:
                error_one_sigma = np.percentile(abs_errors[index_selected], 68.27)
            # Add it to the list
            self._errors_one_sigma.append(self.factor * error_one_sigma)

            # 4 - Build text of water body size (for plotting X labels)
            if (s == nb_slices - 1) and (end_slices_list[s] == "+"):
                if water_body == 'lake size':
                    slice_title = '[>{0}{1}]'.format(int(np.sqrt(cur_slice)), '\u00B2')  # for example: >1000\u00B2
                else:
                    slice_title = '[>{0}]'.format(cur_slice)  # for example: >1000
            else:
                if water_body == 'lake size':
                    slice_title = '[{0}{1}-{2:d}{3}]'.format(int(np.sqrt(cur_slice)), '\u00B2',
                                                             int(np.sqrt(end_slices_list[s])),
                                                             '\u00B2')  # for example: [50\u00B2-100\u00B2]
                else:
                    slice_title = '[{0}-{1}]'.format(cur_slice, end_slices_list[s])  # for example: [50-100]

            # Add slice title to the list
            self._slices_title_list.append(slice_title)

        # 5 - Compute global 1sigma error(s)
        if len(ind_for_global_1) > 0:
            global_1s_error_1 = np.percentile(abs_errors[np.asarray(list(ind_for_global_1))], 68.27)
        else:
            global_1s_error_1 = 0
        if len(ind_for_global_2) > 0:
            global_1s_error_2 = np.percentile(abs_errors[np.asarray(list(ind_for_global_2))], 68.27)
        else:
            global_1s_error_2 = 0
        if len(ind_for_global_3) > 0:
            global_1s_error_3 = np.percentile(abs_errors[np.asarray(list(ind_for_global_3))], 68.27)
        else:
            global_1s_error_3 = 0
        self._errors_one_sigma.insert(0, global_1s_error_1 * self.factor)
        self._errors_one_sigma.insert(1, global_1s_error_2 * self.factor)

        # 6 - Print global mean/median values in log file
        if self._metrix_name.startswith("area_"):
            logger.info("> Compute indexes for global mean/median: ind = ind2")
            ind_for_global = ind_for_global_2
        else:
            logger.info("> Compute indexes for global mean/median: ind = ind2 + ind1")
            ind_for_global = ind_for_global_2.union(ind_for_global_1)
        logger.info(
            "> ind 1 {} ind 2 {} global {}".format(len(ind_for_global_1), len(ind_for_global_2), len(ind_for_global)))

        print("> {0} / {1} / global 1s error slices (1): {2:.03f} {3}".format(self._metrix_name, which_perfos,
                                                                              global_1s_error_1 * self.factor,
                                                                              self.unit))
        print("> {0} / {1} / global 1s error slices (2): {2:.03f} {3}".format(self._metrix_name, which_perfos,
                                                                              global_1s_error_2 * self.factor,
                                                                              self.unit))
        print("> {0} / {1} / global 1s error slices (3): {2:.03f} {3}".format(self._metrix_name, which_perfos,
                                                                              global_1s_error_3 * self.factor,
                                                                              self.unit))

    # ----------------------------------

    def plot_nb_obs_as_a_function_of_water_body_size(self, water_body='lake size'):
        """
        Plot the number of observations by slices of water body size (river witdh or lake size)
        
        Args:
            water_body (str): water body type, among 'river witdh' or 'lake size' (default)
            threshold (float): defined threshold option to ignore abnormal error values
            
        Returns:
            fig: the output figure
        """

        # 1 - Prepare data
        nb_slices = len(self._slices_title_list)  # Nb of slices
        # 1.1 - Compute related array of x indexes
        ind_x = range(nb_slices)
        # 1.2 - Nb of items per slice
        tab_nb_items = np.zeros(nb_slices)
        for cur_slice in ind_x:  # Loop over slices
            tab_nb_items[cur_slice] = len(self._index_sliced[cur_slice][0])
        nb_total = np.sum(tab_nb_items)

        # 2 - Init figure
        fig = plt.figure(figsize=(10, 4))
        plt.ylim(0, np.max(tab_nb_items) * 1.1)

        # 3 - Plot histogram of nb of items vs water body size slice
        barlist = plt.bar(ind_x, tab_nb_items)
        plt.xticks(ind_x, self._slices_title_list)
        # Define a different color for each water bodies size
        if len(self._slices_title_list) == 4:
            colorbar = {0: 'darkviolet', 1: 'violet', 2: 'darkorange', 3: 'gold'}
        else:
            colorbar = {0: 'darkblue', 1: 'darkviolet', 2: 'violet', 3: 'darkorange', 4: 'gold', 5: 'darkkhaki'}
        for ind in ind_x:
            barlist[ind].set_color(colorbar[ind])

        # 4 - Add title and labels
        # 4.1 - Title
        item_name = "observations"
        title = 'Number of matching %s as a function of %s (%d in total)' % (item_name, water_body, nb_total)
        plt.title(title)
        # 4.2 - Label X
        if water_body == 'river width':
            water_body_unit = 'm'
        else:
            water_body_unit = 'm' + '\u00B2'
        xlabel = '{0} ({1})'.format(water_body.capitalize(), water_body_unit)
        plt.xlabel(xlabel)
        # 4.3 - Label Y
        plt.ylabel('Number of %s' % item_name)
        # 4.4 - Nb and associated percentage
        for ind, value in enumerate(tab_nb_items):
            print_str = "{0:.0f}%".format(round(value / nb_total * 100.0))
            plt.text(ind - 0.1, value + 1, print_str, color=colorbar[ind], fontweight='bold')

        return fig

    def plot_along_swath(self, water_body='river width', outlier_lim=None):
        """
        Plot the errors distribution along swath width
        
        Args:
            water_body (str): water body type, among 'lake size' or 'river width'
            outlier_lim (list): define outlier limits. Examples: Area: +/-200 %, WSE: +/-10m
        Returns:
            None
        """
        logger = logging.getLogger(self.__class__.__name__)

        # 1 - Prepare data
        # 1.1 - Errors to plot
        plot_errors = self._errors * self.factor
        # 1.2 - Convert Xtrak_dist to a numpy array
        xrange = np.abs(np.asarray(self._xrange))

        # 2 - Init figure
        fig = plt.figure(figsize=(10, 4))
        plt.subplots_adjust(right=0.75)

        # 2.1 - Define a different color for each water bodies size
        if len(self._slices_title_list) == 4:
            colorbar = {0: 'darkviolet', 1: 'violet', 2: 'darkorange', 3: 'gold'}
        else:
            colorbar = {0: 'darkblue', 1: 'darkviolet', 2: 'violet', 3: 'darkorange', 4: 'gold', 5: 'darkkhaki'}

        # 3 - Plot errors distribution as a function of swath width, distinguish water body size by a different color
        for s in range(len(self._slices_title_list)):  # Loop over slices
            index_selected = self._index_sliced[s]

            outliers = np.zeros(np.shape(xrange[index_selected]), dtype=bool)
            indep = xrange[index_selected].copy()
            dep = plot_errors[index_selected].copy()

            # Clip outliers and plot them with different marker
            if outlier_lim is not None:
                if self._metrix_name == "wse":
                    clipped_value = outlier_lim[1]
                else:
                    clipped_value = outlier_lim[0]
                logger.info('Clip outliers to: +/-{} {}'.format(clipped_value, self.unit))
                outliers = np.logical_or((dep < -clipped_value), (dep > clipped_value))
                dep[dep < -clipped_value] = -clipped_value
                dep[dep > clipped_value] = clipped_value
            not_outliers = np.logical_not(outliers)
            nb_outliers = len(dep[outliers])
            nb_not_outliers = len(dep[not_outliers])
            logger.info('Number of values {} clipped: {} not clipped {}'.format(len(dep), nb_outliers, nb_not_outliers))
            nb_selected = len(index_selected)

            if nb_selected != 0:
                if outlier_lim is None:
                    plt.scatter(indep, dep, c=colorbar[s], label=self._slices_title_list[s])
                else:
                    plt.scatter(indep[not_outliers], dep[not_outliers], c=colorbar[s], label=self._slices_title_list[s])
                    if nb_outliers != 0:
                        plt.scatter(indep[outliers], dep[outliers], c=colorbar[s],
                                    label="(outliers clipped)", marker='+')

        # 3.1 - Dispay legend
        if water_body == 'river width':
            water_body_unit = 'm'
        else:
            water_body_unit = 'm' + '\u00B2'
        legend_title = '{0} ({1})'.format(water_body.capitalize(), water_body_unit)
        plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title=legend_title)

        # 3.2 - Add horizontal line for 0
        if np.min(plot_errors) < 0:
            plt.axhline(0, color="black")

        # 3.3 - Title and labels
        if self._metrix_name == "wse":
            title = 'Error in water surface elevation as a function of swath position'
        else:
            title = 'Relative error in {} surface area as a function of swath position'.format(
                self._metrix_name.split('_')[-1])
        plt.title(title)
        plt.xlabel('Cross-track distance from nadir (km)')
        if self._metrix_name == "wse":
            ylabel = 'WSE error ({0})'.format(self.unit)
        else:
            ylabel = 'Relative area error ({0})'.format(self.unit)
        plt.ylabel(ylabel)

        return fig

    def plot_along_swath_popular(self, threshold=None):
        """
        Plot the errors distribution along swath width
        Args:
            threshold (float): defined threshold option to ignore abnormal error values
        Returns:
            None
        """

        logger = logging.getLogger(self.__class__.__name__)

        # 1 - Prepare data
        # 1.1 - Copy of tab of errors
        plot_errors = self._errors.copy()
        # 1.2 - Convert Xtrak_dist to a numpy array
        xrange = np.abs(np.asarray(self._xrange))
        # 1.3 - Retrieve max error
        max_error = np.max(plot_errors) * self.factor
        logger.info('plot_errors MAX: {}'.format(max_error))

        ##xrange = np.round(xrange/5.)*5.
        # xrange = np.round((xrange+2.55)/5.)*5.-2.5 ## POUR L'ARTICLE

        # 2 - Build DataFrame wrt to threshold or not
        if threshold is not None:
            index_threshold = np.where(plot_errors < threshold)
            df = pd.DataFrame({'y': plot_errors[index_threshold] * self.factor, 'x': xrange[index_threshold]})
        else:
            df = pd.DataFrame({'y': plot_errors * self.factor, 'x': xrange})
        # indexes = df.index
        df = df.groupby('x').median()

        # 3 - Init figure
        fig = plt.figure(figsize=(10, 4))
        # plt.bar(indexes, df['y'])
        plt.bar(df.index, df['y'])
        # plt.bar(df.index, df['y'], width=4.5) ## POUR L'ARTICLE
        # 4 - Plot measured and required 1sigma error
        # 4.1 - Display measured relative error one sigma
        # plt.axhline(self._errors_one_sigma[0], color="blue", linestyle="--",
        #            label="Measured $1\sigma$ (%s)" % self.print_metrix_value(self._errors_one_sigma[0]))
        # 4.2 - Display requirement value
        # plt.axhline(self.err_req, color="red",
        #            label="Requirement (%s)" % self.print_metrix_value(self.err_req))
        # 4.3 - Add legend
        # plt.legend()
        # 4.4 - Add horizontal line for 0
        if np.min(plot_errors) < 0:
            plt.axhline(0, color="black")

        # 4 - Add title and labels
        # 4.2 - Label X
        plt.xlabel('Cross-track distance from nadir (km)')
        # 4.3 - Label Y
        if self._metrix_name == "wse":
            plt.ylabel('Median WSE error (m)')
            title = 'Median error in water surface elevation as a function of swath position'
        else:
            plt.ylabel('Median relative area error (%)')
            title = 'Median relative error in {} surface area as a function of swath position'.format(
                self._metrix_name.split('_')[-1])
        # plt.title(title)

        return fig

    def plot_one_sigma_error_distribution(self, water_body='river width', which_perfos='req'):
        """
        Plot the measurement error one sigma by slices of water body size

        Args:
            water_body (str): water body type, among 'lake size' or 'river width'
            which_perfos (str): 'req' or 'goal'
        Returns:
            None
        """
        logger = logging.getLogger(self.__class__.__name__)
        # 1 - Prepare data
        # 1.1 - Convert 1sigma error list to numpy array
        one_sigma_error_tab = np.asarray(self._errors_one_sigma[2:])
        # 1.2 - Compute related array of x indexes
        ind_x = np.arange(len(self._errors_one_sigma[2:]))

        # 2 - Init figure
        fig = plt.figure(figsize=(10, 4))
        # 2.1 - Plot histogram of 1sigma error vs water body size slice
        barlist = plt.bar(ind_x, one_sigma_error_tab)
        plt.xticks(ind_x, self._slices_title_list)
        # 2.2 - Define a different color for each water bodies size
        if len(self._slices_title_list) == 4:
            colorbar = {0: 'darkviolet', 1: 'violet', 2: 'darkorange', 3: 'gold'}
        else:
            colorbar = {0: 'darkblue', 1: 'darkviolet', 2: 'violet', 3: 'darkorange', 4: 'gold', 5: 'darkkhaki'}
        for ind in ind_x:
            barlist[ind].set_color(colorbar[ind])

        # 3 - Plot measured and required 1sigma error
        if which_perfos == 'goal':
            x_min = 1 / len(ind_x) * 1.1
            x_th = 1 / len(ind_x) * 3
        else:
            x_min = 0
            x_th = 1 / len(ind_x) * 2
        if self._metrix_name.startswith("area_"):
            # 3.1 - Display measured relative error one sigma:
            # only measured > self._req_th_legend (100m for rivers, 250²m² for lakes)
            plt.axhline(self._errors_one_sigma[1], xmin=x_min, xmax=1, color="blue", linestyle="--",
                        label="Measured > %s (%.3f)" % (
                            self._req_th_legend, self._errors_one_sigma[1]))
            # 3.2 - Display requirement value
            plt.axhline(y=self.err_req, xmin=x_min, xmax=1, color="red",
                        label="Requirement > %s (%s)" % (
                            self._req_th_legend, self.print_metrix_value(self.err_req)))
            if which_perfos == 'goal':
                # 3.3 - Display goal value
                label = "Goal < %s (%s)" % (
                    self._req_th_legend, self.print_metrix_value(self.err_goal))
                plt.axhline(y=self.err_goal, xmin=0, xmax=x_min, color="red", linestyle=":", label=label)
        else:
            # 3.1 - Display measured value for lakes surface < 1km²
            label = "Measured < %s (%s)" % (self._req_th_legend, self.print_metrix_value(self._errors_one_sigma[0]))
            plt.axhline(self._errors_one_sigma[0], xmin=x_min, xmax=x_th, color="blue", linestyle=":", label=label)
            # 3.2 - Display requirement value for lakes surface < 1km²
            label = "Requirement < %s (%s)" % (self._req_th_legend, self.print_metrix_value(self.err_goal))
            plt.axhline(y=self.err_goal, xmin=x_min, xmax=x_th, color="red", linestyle=":", label=label)
            # 3.3 - Display measured value for lakes surface > 1km²
            label = "Measured > %s (%s)" % (self._req_th_legend, self.print_metrix_value(self._errors_one_sigma[1]))
            plt.axhline(self._errors_one_sigma[1], xmin=x_th, xmax=1, color="blue", linestyle="--", label=label)
            # 3.4 - Display requirement value for lakes surface > 1km²
            label = "Requirement > %s (%s)" % (self._req_th_legend, self.print_metrix_value(self.err_req))
            plt.axhline(y=self.err_req, xmin=x_th, xmax=1, color="red", label=label)
        logging.info("Error 1 sigma for %s for %s : %f %f" % (
        self._metrix_name, which_perfos, self._errors_one_sigma[0], self._errors_one_sigma[1]))
        # 3.4 - Add legend
        plt.legend()

        # 4 - Add title and labels
        # 4.1 - Title
        if self._metrix_name == "wse":
            title = '|Error| in %s (1$\sigma$) as a function of %s' % (self.metrix_long_name, water_body)
        else:
            title = '|%s error| in %s %s (1$\sigma$) as a function of %s' % (
                self.error_type, self._metrix_name.split('_')[-1], self.metrix_long_name, water_body)
        plt.title(title)
        # 4.2 - Label X
        if water_body == 'river width':
            water_body_unit = 'm'
        else:
            water_body_unit = 'm' + '\u00B2'
        xlabel = '{0} ({1})'.format(water_body.capitalize(), water_body_unit)
        plt.xlabel(xlabel)
        # 4.3 - Label Y
        if self._metrix_name == "wse":
            ylabel = '$1\sigma$ |WSE error| ({0})'.format(self.unit)
        else:
            ylabel = '$1\sigma$ |relative area error| ({0})'.format(self.unit)
        plt.ylabel(ylabel)

        return fig
