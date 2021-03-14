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

	const auto results = Simulate(params);

	const auto [uniform_mean, uniform_variance] = uniform_mean_and_variance(params.a, params.b);
	const auto [weibull_mean, weibull_variance] = weibull_mean_and_variance(params.k, params.lambda);
	const auto rho =  weibull_mean / uniform_mean;

	ui->muUniformLineEdit->setText(QString::number(uniform_mean));
	ui->sigma2UniformLineEdit->setText(QString::number(uniform_variance));
	ui->muWeibullLineEdit->setText(QString::number(weibull_mean));
	ui->sigma2WeibullLineEdit->setText(QString::number(weibull_variance));

	ui->estimatedSystemLoadLineEdit->setText(QString::number(rho));
	ui->actualSystemLoadLineEdit->setText(QString::number(results.load));
}
