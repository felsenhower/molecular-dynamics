/******************************************************************************/
/*                                                                            */
/*                             Molecular dynamics                             */
/*                             ==================                             */
/*                                                                            */
/* Copyright (c) 2017 Ruben Felgenhauer, Leonhard Reichenbach                 */
/*                                                                            */
/* This file is released under the MIT license, but is made to be used with   */
/* the Qt Framework, which is released under the GNU Lesser General Public    */
/* License (LGPL) and partly under the GNU General Public License (GPL).      */
/*                                                                            */
/******************************************************************************/

#include "mainwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
  QApplication a(argc, argv);

  uint64_t number_of_threads = 1;
  uint64_t N = 25;
  double L = 5.0;
  double timestep = 0.005;
  double max_deviation = 0.0001;
  double max_v0 = 0.000;
  bool periodic_boundaries = true;
  double dist_threshold = 3.0;

  SimParams sim_params = init_params(number_of_threads, N, L, timestep, max_deviation,
                                     max_v0, periodic_boundaries, dist_threshold, 0);

  DisplayOptions display_options = {.frameskip = 0,
                                    .dist_threshold_indicator_enabled = false,
                                    .dist_threshold_color = Qt::blue,
                                    .sigma_indicator_enabled = false,
                                    .sigma_color = Qt::red,
                                    .particle_size = 10,
                                    .update_interval = 1.0};

  MainWindow w(NULL, sim_params, display_options);
  w.show();

  return a.exec();
}
