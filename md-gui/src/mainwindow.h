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

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include <QtGui>
#include <QtCore>
#include <QColorDialog>
#include <QFrame>
#include <QMessageBox>

extern "C" {
    #include "../../md-core/md-core.h"
}

typedef struct {
  int frameskip;
  bool dist_threshold_indicator_enabled;
  QColor dist_threshold_color;
  bool sigma_indicator_enabled;
  QColor sigma_color;
  int particle_size;
  double update_interval;
} DisplayOptions;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  explicit MainWindow(QWidget *parent = 0);
  MainWindow(QWidget *parent, SimParams sim_params, DisplayOptions display_options);
  ~MainWindow();
  void rescaleUi();
  void print();
  void set_ui_bindings();
  SimParams sim_params;
  SimData sim_data;
  struct timeval before_step;
  struct timeval after_step;
  bool do_step = false;
  double update_counter = 0.0;
  void step_ui();
  DisplayOptions display_options;
  QRect getRect(QFrame *frame);
  void drawSimulation();
  void drawGraph(QFrame *canvas, std::vector<double> values, const QString &label);
  std::vector<double> energies;
  std::vector<double> temperatures;

protected:
  //void resizeEvent(QResizeEvent *);
  void paintEvent(QPaintEvent *);

private slots:
  void increment_spinbox();
  void decrement_spinbox();
  void next_sqr_spinbox();
  void prev_sqr_spinbox();
  void focus_spinbox();
  void set_frm_color();
  void on_btn_sim_params_apply_clicked();
  void on_btn_sim_step_clicked();
  void on_btn_sim_start_clicked();
  void on_btn_sim_stop_clicked();
  void on_btn_display_apply_clicked();
  void on_btn_heat_apply_clicked();

  void on_btn_heat_abort_clicked();
  
  void on_btn_help_clicked();
  
private:
  Ui::MainWindow *ui;

};

#endif // MAINWINDOW_H
