#include "mainwindow.hpp"
#include "./ui_mainwindow.h"

#include "random.hpp"
#include "simulate.hpp"

MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
	, ui(new Ui::MainWindow)
{
	ui->setupUi(this);
}

MainWindow::~MainWindow()
{
	delete ui;
}

void MainWindow::on_simulatePushButton_clicked()
{
	const SimulateParams params = {
		.a = ui->aDoubleSpinBox->value(),
		.b = ui->bDoubleSpinBox->value(),

		.k = ui->kDoubleSpinBox->value(),
		.lambda = ui->lambdaDoubleSpinBox->value(),

		.t = ui->tDoubleSpinBox->value(),
	};

	const auto [uniform_mean, uniform_variance] = uniform_mean_and_variance(params.a, params.b);
	const auto [weibull_mean, weibull_variance] = weibull_mean_and_variance(params.k, params.lambda);
	const auto rho =  weibull_mean / uniform_mean;

	const auto result = Simulate(params);

	ui->uniformMeanLineEdit->setText(QString::number(uniform_mean));
	ui->uniformVarianceLineEdit->setText(QString::number(uniform_variance));

	ui->weibullMeanLineEdit->setText(QString::number(weibull_mean));
	ui->weibullVarianceLineEdit->setText(QString::number(weibull_variance));

	ui->inputIntensityLineEdit->setText(QString::number(1.0 / uniform_mean));
	ui->outputIntensityLineEdit->setText(QString::number(1.0 / weibull_mean));

	ui->estimatedSystemLoadLineEdit->setText(QString::number(rho));
	ui->actualSystemLoadLineEdit->setText(QString::number(result.load));
}

void MainWindow::on_simulate2PushButton_clicked()
{
	const auto input_intensity = ui->inputIntensity2DoubleSpinBox->value();
	const auto uniform_std = ui->uniformStd2DoubleSpinBox->value();
	const auto output_intensity = ui->outputIntensity2DoubleSpinBox->value();
	const auto weibull_std = ui->weibullStd2DoubleSpinBox->value();
	const auto t = ui->t2DoubleSpinBox->value();

	const auto uniform_mean = 1 / input_intensity;
	const auto weibull_mean = 1 / output_intensity;

	const auto [a, b] = uniform_parameters_from_mean_and_std(uniform_mean, uniform_std);
	const auto [k, lambda] = weibull_parameters_from_mean_and_std(weibull_mean, weibull_std);

	const auto rho =  weibull_mean / uniform_mean;

	const SimulateParams params = {
		.a = a,
		.b = b,

		.k = k,
		.lambda = lambda,

		.t = t,
	};

	const auto result = Simulate(params);

	ui->uniformMean2LineEdit->setText(QString::number(uniform_mean));
	ui->weibullMean2LineEdit->setText(QString::number(weibull_mean));

	ui->a2LineEdit->setText(QString::number(a));
	ui->b2LineEdit->setText(QString::number(b));

	ui->k2LineEdit->setText(QString::number(k));
	ui->lambda2LineEdit->setText(QString::number(lambda));

	ui->estimatedSystemLoadLineEdit->setText(QString::number(rho));
	ui->actualSystemLoadLineEdit->setText(QString::number(result.load));
}
