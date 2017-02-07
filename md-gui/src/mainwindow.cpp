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
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::MainWindow)
{
  ui->setupUi(this);
  rescaleUi();
}

MainWindow::MainWindow(QWidget *parent, SimParams sim_params, DisplayOptions display_options) :
  QMainWindow(parent),
  ui(new Ui::MainWindow)
{
  this->sim_params = sim_params;
  this->display_options = display_options;
  
  seed();

  this->sim_data = init_data(&sim_params);

  gettimeofday(&(this->before_step), NULL);
  gettimeofday(&(this->after_step), NULL);

  ui->setupUi(this);

  this->set_ui_bindings();

  ui->sbox_num_threads->setValue(this->sim_params.number_of_threads);
  ui->sbox_num_particles->setValue(this->sim_params.N);
  ui->sbox_lattice_len->setValue(this->sim_params.L);
  ui->sbox_timestep->setValue(this->sim_params.timestep);
  ui->sbox_max_deviation->setValue(this->sim_params.max_deviation);
  ui->sbox_max_v0->setValue(this->sim_params.max_v0);
  ui->cbox_periodic_boundaries->setChecked(this->sim_params.periodic_boundaries);
  ui->sbox_dist_threshold->setValue(this->sim_params.dist_threshold);

  ui->sbox_frameskip->setValue(this->display_options.frameskip);
  ui->cbox_dist_threshold->setChecked(
        this->display_options.dist_threshold_indicator_enabled);
  ui->cbox_sigma->setChecked(this->display_options.sigma_indicator_enabled);

  QFrame *frame; QPalette palette;
  frame = ui->frm_dist_threshold_col;
  palette = frame->palette();
  palette.setColor(QPalette::Window, this->display_options.dist_threshold_color);
  frame->setPalette(palette);
  frame = ui->frm_sigma_col;
  palette = frame->palette();
  palette.setColor(QPalette::Window, this->display_options.sigma_color);
  frame->setPalette(palette);

  ui->sbox_particle_size->setValue(this->display_options.particle_size);
  ui->sbox_stat_update->setValue(this->display_options.update_interval);

  //this->print();
  
}

void MainWindow::set_ui_bindings() {
  std::vector<QPushButton*> plus_buttons = {this->ui->btn_num_threads_p,
                                            this->ui->btn_num_particles_p,
                                            this->ui->btn_lattice_len_p,
                                            this->ui->btn_timestep_p,
                                            this->ui->btn_max_deviation_p,
                                            this->ui->btn_max_v0_p,
                                            this->ui->btn_dist_threshold_p,
                                            this->ui->btn_frameskip_p,
                                            this->ui->btn_particle_size_p,
                                            this->ui->btn_stat_update_p,
                                            this->ui->btn_target_temp_p,
                                            this->ui->btn_heat_time_p/*,
                                            this->ui->btn_attack_coeff_p*/};
  std::vector<QPushButton*> minus_buttons = {this->ui->btn_num_threads_m,
                                            this->ui->btn_num_particles_m,
                                            this->ui->btn_lattice_len_m,
                                            this->ui->btn_timestep_m,
                                            this->ui->btn_max_deviation_m,
                                            this->ui->btn_max_v0_m,
                                            this->ui->btn_dist_threshold_m,
                                            this->ui->btn_frameskip_m,
                                            this->ui->btn_particle_size_m,
                                            this->ui->btn_stat_update_m,
                                            this->ui->btn_target_temp_m,
                                            this->ui->btn_heat_time_m/*,
                                            this->ui->btn_attack_coeff_m*/};
  for(auto it=plus_buttons.begin() ; it < plus_buttons.end(); it++ ) {
    connect(*it, SIGNAL(pressed()), this, SLOT(increment_spinbox()));
    connect(*it, SIGNAL(released()), this, SLOT(focus_spinbox()));
  }
  for(auto it=minus_buttons.begin() ; it < minus_buttons.end(); it++ ) {
    connect(*it, SIGNAL(pressed()), this, SLOT(decrement_spinbox()));
    connect(*it, SIGNAL(released()), this, SLOT(focus_spinbox()));
  }

  connect(this->ui->btn_num_particles_mm, SIGNAL(pressed()),
          this, SLOT(prev_sqr_spinbox()));
  connect(this->ui->btn_num_particles_mm, SIGNAL(released()),
          this, SLOT(focus_spinbox()));
  connect(this->ui->btn_num_particles_pp, SIGNAL(pressed()),
          this, SLOT(next_sqr_spinbox()));
  connect(this->ui->btn_num_particles_pp, SIGNAL(released()),
          this, SLOT(focus_spinbox()));

  connect(this->ui->btn_dist_threshold_col, SIGNAL(clicked()),
          this, SLOT(set_frm_color()));
  connect(this->ui->btn_sigma_col, SIGNAL(clicked()),
          this, SLOT(set_frm_color()));
}

MainWindow::~MainWindow()
{
  this->do_step = false;
  finalize(&(this->sim_data), &(this->sim_params));
  delete ui;
}

void MainWindow::print() {
  for (size_t i = 0; i < sim_params.N; i++) {
    printf(/*foo, */"i=%0ld:\t%.5g\t%.5g\n", i,
            this->sim_data.molecules[i].curr_pos.x,
            this->sim_data.molecules[i].curr_pos.y);
  }
}

void MainWindow::rescaleUi() {
  int height = this->ui->cont_lvl2_1->height();
  int width = this->ui->cont_lvl2_1->width();
  if (height > width) {
    this->ui->spacer_draw_area_left->setMaximumWidth(0);
    this->ui->spacer_draw_area_right->setMaximumWidth(0);
    this->ui->spacer_draw_area_top->setMaximumHeight((height-width)/2);
    this->ui->spacer_draw_area_bottom->setMaximumHeight((height-width)/2);
  } else {
    this->ui->spacer_draw_area_left->setMaximumWidth((width-height)/2);
    this->ui->spacer_draw_area_right->setMaximumWidth((width-height)/2);
    this->ui->spacer_draw_area_top->setMaximumHeight(0);
    this->ui->spacer_draw_area_bottom->setMaximumHeight(0);
  }
}

void MainWindow::step_ui() {
  struct timeval temp;
  struct timeval start_to_start;
  struct timeval start_to_end;

  gettimeofday(&temp, NULL);
  timersub(&temp, &(this->before_step), &start_to_start);
  this->before_step = temp;

  for (int i = 0; i <= this->display_options.frameskip; i++) {
    step(&sim_data, &sim_params);
  }

  gettimeofday(&(this->after_step), NULL);

  timersub(&after_step, &before_step, &start_to_end);

  //this->print();

  this->repaint();

  double start_to_start_seconds = (double)start_to_start.tv_sec
      + (((double)start_to_start.tv_usec) / 1000000);

  double start_to_end_seconds = (double)start_to_end.tv_sec
      + (((double)start_to_end.tv_usec) / 1000000);

  this->update_counter += start_to_start_seconds;
  if (this->update_counter > this->display_options.update_interval) {
    this->update_counter = 0.0;

    this->energies.push_back(this->sim_data.energy);
    if (energies.size() > 1024) {
      energies.erase(energies.begin());
    }
    
    this->temperatures.push_back(this->sim_data.temperature);
    if (temperatures.size() > 1024) {
      temperatures.erase(temperatures.begin());
    }

    this->ui->lbl_calc_time->setText("Calculation time: "
        + QString::number(start_to_end_seconds*1000000, 'f', 0)
        + " µs ");

    this->ui->lbl_frame_time->setText("Frame time: "
        + QString::number(start_to_start_seconds*1000000, 'f', 0)
        + " µs");

    this->ui->lbl_frame_rate->setText("Refresh rate: "
        + QString::number(1/start_to_start_seconds, 'f', 1)
        + " fps");
  }
  
  this->ui->btn_heat_apply->setEnabled(!this->sim_data.heat_bath_enabled);
  this->ui->btn_heat_abort->setEnabled(this->sim_data.heat_bath_enabled);
  
}

QRect MainWindow::getRect(QFrame *frame) {
  int width = frame->contentsRect().width();
  int height = frame->contentsRect().height();
  int x = 0;
  int y = 0;
  QWidget *widget = frame;
  while (widget != (QWidget*)(this)) {
    x += widget->x();
    y += widget->y();
    widget = widget->parentWidget();
  }
  return QRect(x,y,width,height);
}

void MainWindow::drawSimulation() {
  QFrame *canvas = this->ui->draw_area;
  QRect canvasRect = getRect(canvas);

  QPainter painter(this);
  painter.setClipping(true);
  painter.setClipRect(canvasRect);
  painter.fillRect(painter.clipBoundingRect(), Qt::white);

  QPen pointpen(Qt::black);
  pointpen.setWidth(this->display_options.particle_size);
  pointpen.setCapStyle(Qt::RoundCap);

  QPen circlepen(this->display_options.dist_threshold_color);
  circlepen.setCapStyle(Qt::RoundCap);
  circlepen.setWidth(2);

  QPen smallcirclepen(this->display_options.sigma_color);
  circlepen.setCapStyle(Qt::RoundCap);
  circlepen.setWidth(2);

  for (size_t i = 0; i < sim_params.N; i++) {
    QPoint p;

    int x = std::ceil((double)(canvasRect.left()) + (double)(canvasRect.width())
                      * (this->sim_data.molecules[i].curr_pos.x)
                      / ((double)(this->sim_params.L)));

    int y = std::ceil((double)(canvasRect.top()) + (double)(canvasRect.height())
                      * (this->sim_data.molecules[i].curr_pos.y)
                      / ((double)(this->sim_params.L)));

    p.setX(x);
    p.setY(y);

    painter.setPen(pointpen);
    painter.drawPoint(p);

    if (display_options.sigma_indicator_enabled) {
      painter.setPen(smallcirclepen);
      int r = std::ceil(((double)(canvasRect.width())) /
                        ((double)(this->sim_params.L)));
      painter.drawEllipse(p, r, r);
    }

    if (display_options.dist_threshold_indicator_enabled) {
      painter.setPen(circlepen);
      int r = std::ceil(((double)(canvasRect.width())) *
                        this->sim_params.dist_threshold /
                        ((double)(this->sim_params.L)));
      painter.drawEllipse(p, r, r);
    }
  }
}

void MainWindow::drawGraph(QFrame *canvas, std::vector<double> values,
                           const QString &label) {
  QRect canvasRect = getRect(canvas);
  QPainter painter(this);
  painter.setClipping(true);
  painter.setClipRect(canvasRect);
  painter.fillRect(painter.clipBoundingRect(), Qt::white);
  
  if (!values.empty()) {
    QPen linepen(Qt::black);
    QPen pointpen(Qt::black);
    pointpen.setCapStyle(Qt::RoundCap);
    pointpen.setWidth(5);

    double max_value = *(std::max_element(values.begin(), values.end()));
    double min_value = *(std::min_element(values.begin(), values.end()));

    int right = canvasRect.right() - 215;
    int bottom = canvasRect.top() + canvasRect.height() - 15;
    int height = canvasRect.height() - 30;
    int width = canvasRect.width();
    int step = 10;

    int i = 0;
    QPoint old_point;
    QPoint new_point;
    for (auto it = values.rbegin();
         it < values.rend() && i * step < width;
         it++, i++) {
      double value = *it;
      int x = right - i * step;
      int y = ceil((double)bottom - (value - min_value) *
                   (double)height / (max_value - min_value));
      if (!i) {
        old_point.setX(x);
        old_point.setY(y);
        painter.setPen(pointpen);
        painter.drawPoint(old_point);
        painter.setPen(linepen);
        painter.drawText(old_point.x()+15,old_point.y(),
                         label + ": " + QString::number(value, 'g', 4));
      } else {
        new_point.setX(x);
        new_point.setY(y);
        painter.drawLine(old_point, new_point);
        old_point.setX(new_point.x());
        old_point.setY(new_point.y());
      }
    }
  }
}

void MainWindow::paintEvent(QPaintEvent *) {
  rescaleUi();

  drawSimulation();

  drawGraph(this->ui->frm_energy_graph, this->energies, "Energy");
  
  drawGraph(this->ui->frm_temp_graph, this->temperatures, "Temperature");
}

void MainWindow::decrement_spinbox() {
  QAbstractSpinBox *spinbox = ((QPushButton *)(sender()))->parentWidget()
                    ->findChild<QAbstractSpinBox*>();
  QSpinBox *_spinbox = qobject_cast<QSpinBox*>(spinbox);
  if (_spinbox) {
    _spinbox->setValue(_spinbox->value() - _spinbox->singleStep());
  } else {
    QDoubleSpinBox *_spinbox = qobject_cast<QDoubleSpinBox*>(spinbox);
    _spinbox->setValue(_spinbox->value() - _spinbox->singleStep());
  }
}

void MainWindow::increment_spinbox() {
  QAbstractSpinBox *spinbox = ((QPushButton *)(sender()))->parentWidget()
                    ->findChild<QAbstractSpinBox*>();
  QSpinBox *_spinbox = qobject_cast<QSpinBox*>(spinbox);
  if (_spinbox) {
    _spinbox->setValue(_spinbox->value() + _spinbox->singleStep());
  } else {
    QDoubleSpinBox *_spinbox = qobject_cast<QDoubleSpinBox*>(spinbox);
    _spinbox->setValue(_spinbox->value() + _spinbox->singleStep());
  }
}

void MainWindow::prev_sqr_spinbox() {
  QSpinBox *spinbox = ((QPushButton *)(sender()))->parentWidget()
                    ->findChild<QSpinBox*>();
  int temp = (ceil(sqrt((double)spinbox->value()))-1);
  spinbox->setValue(temp*temp);
}

void MainWindow::next_sqr_spinbox() {
  QSpinBox *spinbox = ((QPushButton *)(sender()))->parentWidget()
                    ->findChild<QSpinBox*>();
  int temp = (floor(sqrt((double)spinbox->value()))+1);
  spinbox->setValue(temp*temp);
}

void MainWindow::focus_spinbox() {
  QPushButton *_sender = (QPushButton *)sender();
  if (!_sender->isDown()) {
  QAbstractSpinBox *spinbox = _sender->parentWidget()
                    ->findChild<QAbstractSpinBox*>();
  spinbox->setFocus();
  }
}

void MainWindow::set_frm_color() {
  QFrame *frame = (QFrame *)(((QPushButton*)(sender()))->parentWidget());
  QPalette palette = frame->palette();
  QColor color = QColorDialog::getColor(palette.window().color(), this);
  if(color.isValid()) {
    palette.setColor(QPalette::Window, color);
  }
  frame->setPalette(palette);
}

void MainWindow::on_btn_sim_params_apply_clicked()
{
  finalize(&(this->sim_data), &(this->sim_params));
  this->energies.clear();
  this->temperatures.clear();
  this->sim_params = init_params(ui->sbox_num_threads->value(),
                                 ui->sbox_num_particles->value(),
                                 ui->sbox_lattice_len->value(),
                                 ui->sbox_timestep->value(),
                                 ui->sbox_max_deviation->value(),
                                 ui->sbox_max_v0->value(),
                                 ui->cbox_periodic_boundaries->isChecked(),
                                 ui->sbox_dist_threshold->value(),
                                 0);
  this->sim_data = init_data(&sim_params);
  //this->print();
  this->repaint();
}

void MainWindow::on_btn_sim_step_clicked()
{
  this->step_ui();
}

void MainWindow::on_btn_sim_start_clicked()
{
  this->do_step = true;
  this->ui->btn_sim_start->setEnabled(false);
  this->ui->btn_sim_stop->setEnabled(true);
  this->ui->btn_sim_step->setEnabled(false);
  while (this->do_step && this->isVisible()) {
    this->step_ui();
    QApplication::processEvents();
  }
}

void MainWindow::on_btn_sim_stop_clicked()
{
  this->ui->btn_sim_start->setEnabled(true);
  this->ui->btn_sim_stop->setEnabled(false);
  this->ui->btn_sim_step->setEnabled(true);
  this->do_step = false;
}

void MainWindow::on_btn_display_apply_clicked()
{
  this->display_options.frameskip = this->ui->sbox_frameskip->value();
  this->display_options.dist_threshold_indicator_enabled = this->ui->cbox_dist_threshold->isChecked();
  this->display_options.dist_threshold_color = this->ui->frm_dist_threshold_col->palette().window().color();
  this->display_options.sigma_indicator_enabled = this->ui->cbox_sigma->isChecked();
  this->display_options.sigma_color = this->ui->frm_sigma_col->palette().window().color();
  this->display_options.particle_size = this->ui->sbox_particle_size->value();
  this->display_options.update_interval = this->ui->sbox_stat_update->value();
  this->repaint();
}

void MainWindow::on_btn_heat_apply_clicked()
{
  double target_temperature = this->ui->sbox_target_temp->value();
  uint64_t steps = 100 * this->ui->sbox_heat_time->value();
  //double attack_coefficient = this->ui->sbox_attack_coeff->value();

  enable_heat_bath(&(this->sim_data), target_temperature, steps);
  
  this->ui->btn_heat_apply->setEnabled(false);
  this->ui->btn_heat_abort->setEnabled(true);
}

void MainWindow::on_btn_heat_abort_clicked()
{
  this->sim_data.heat_bath_enabled = false;  
  this->ui->btn_heat_apply->setEnabled(true);
}

void MainWindow::on_btn_help_clicked()
{
  QMessageBox msg;
  msg.setInformativeText("Copyright (c) 2017 Ruben Felgenhauer, Leonhard "
                         "Reichenbach" "\n\n" "This application is released "
                         "under the MIT license, but is made with the Qt "
                         "Framework, which is released under the GNU Lesser "
                         "General Public License (LGPL) and partly under the " 
                         "GNU General Public License (GPL).");
  msg.setText("Molecular dynamics");
  msg.setStandardButtons(QMessageBox::Close);
  msg.exec();
}
