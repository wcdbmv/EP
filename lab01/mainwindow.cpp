#include "mainwindow.hpp"
#include "./ui_mainwindow.h"

#include "random.hpp"
#include "simulate.hpp"

void setTitleOfPlot(QCustomPlot* customPlot, const QString& title) {
	customPlot->addGraph();
	customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

	customPlot->plotLayout()->insertRow(0);
	customPlot->plotLayout()->addElement(0, 0, new QCPTextElement(customPlot, title, QFont("sans", 12, QFont::Bold)));

	customPlot->xAxis->setLabel("ρ");
	customPlot->yAxis->setLabel("Среднее время ожидания");

	customPlot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignLeft | Qt::AlignTop);
}

MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
	, ui(new Ui::MainWindow)
{
	ui->setupUi(this);

	setTitleOfPlot(ui->customPlot, "Cреднее время ожидания (при lambda = 5)");
	setTitleOfPlot(ui->customPlot2, "Среднее время ожидания (при mu = 5)");
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

	plot();
}

void MainWindow::on_simulate2PushButton_clicked()
{
	const auto input_intensity = ui->inputIntensity2DoubleSpinBox->value();
	const auto uniform_std = ui->uniformStd2DoubleSpinBox->value();
	const auto output_intensity = ui->outputIntensity2DoubleSpinBox->value();
	const auto weibull_std = ui->weibullStd2DoubleSpinBox->value();
	const auto t = ui->t2DoubleSpinBox->value();

	const auto uniform_mean = 1.0 / input_intensity;
	const auto weibull_mean = 1.0 / output_intensity;

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

	plot();
}

void MainWindow::plot() {
	constexpr double input_intensity_1 = 5.0;
	constexpr double uniform_mean_1 = 1.0 / input_intensity_1;

	constexpr double output_intensity_2 = 5.0;
	constexpr double weibull_mean_2 = 1.0 / output_intensity_2;

	constexpr double rho_min = 0.1;
	constexpr double rho_max = 0.99;
	constexpr double drho = 0.01;

	constexpr double uniform_std = 0.05;
	constexpr double weibull_std = 0.05;

	constexpr double t = 1000.0;

	QVector<double> x;
	QVector<double> y1;
	QVector<double> y2;

	for (double rho = rho_min; rho <= rho_max; rho += drho) {
		const double output_intensity_1 = input_intensity_1 / rho;
		const double weibull_mean_1 = 1.0 / output_intensity_1;

		const auto [a_1, b_1] = uniform_parameters_from_mean_and_std(uniform_mean_1, uniform_std);
		const auto [k_1, lambda_1] = weibull_parameters_from_mean_and_std(weibull_mean_1, weibull_std);

		const SimulateParams params_1 = {
			.a = a_1,
			.b = b_1,

			.k = k_1,
			.lambda = lambda_1,

			.t = t,
		};

		const auto result_1 = Simulate(params_1);

		const double input_intensity_2 = rho * output_intensity_2;
		const double uniform_mean_2 = 1.0 / input_intensity_2;

		const auto [a_2, b_2] = uniform_parameters_from_mean_and_std(uniform_mean_2, uniform_std);
		const auto [k_2, lambda_2] = weibull_parameters_from_mean_and_std(weibull_mean_2, weibull_std);

		const SimulateParams params_2 = {
			.a = a_2,
			.b = b_2,

			.k = k_2,
			.lambda = lambda_2,

			.t = t,
		};

		const auto result_2 = Simulate(params_2);

		x.push_back(rho);
		y1.push_back(result_1.average_waiting);
		y2.push_back(result_2.average_waiting);
	}

	ui->customPlot->graph(0)->setData(x, y1);
	ui->customPlot2->graph(0)->setData(x, y2);

	ui->customPlot->rescaleAxes();
	ui->customPlot2->rescaleAxes();

	ui->customPlot->replot();
	ui->customPlot2->replot();
}
