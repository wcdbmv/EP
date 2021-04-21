#pragma once

#include <QMainWindow>

#include "experiment.hpp"


QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow(QWidget *parent = nullptr);
	~MainWindow();

private slots:
	void on_calculatePushButton_clicked();

private:
	Ui::MainWindow *ui;

	static constexpr FfeParameters FFE_PARAMS{
		.lambda1       = {.min = 1.0,  .max = 1.5},
		.lambda2       = {.min = 1.0,  .max = 1.5},
		.mu1           = {.min = 5.0,  .max = 6.0},
		.mu2           = {.min = 5.0,  .max = 6.0},
		.sigma_lambda1 = {.min = 0.01, .max = 0.05},
		.sigma_lambda2 = {.min = 0.01, .max = 0.05},
		.times = 100,
	};

	FfeResult full_result;
	FfeResult frac_result;
};
