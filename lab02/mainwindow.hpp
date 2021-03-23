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

	static constexpr double LAMBDA_MIN = 0.5;
	static constexpr double LAMBDA_MAX = 3.0;
	static constexpr double SIGMA_LAMBDA_MIN = 0.01;
	static constexpr double SIGMA_LAMBDA_MAX = 0.10;
	static constexpr double MU_MIN = 5.0;
	static constexpr double MU_MAX = 6.0;

	static constexpr FfeParameters PARAMS{
		.lambda_min = LAMBDA_MIN,
		.lambda_max = LAMBDA_MAX,
		.sigma_lambda_min = SIGMA_LAMBDA_MIN,
		.sigma_lambda_max = SIGMA_LAMBDA_MAX,
		.mu_min = MU_MIN,
		.mu_max = MU_MAX,
		.times = 5,
	};

	FfeResult result;
};
