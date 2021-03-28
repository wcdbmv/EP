#include "mainwindow.hpp"
#include "./ui_mainwindow.h"

#include "random.hpp"


MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
	, ui(new Ui::MainWindow)
{
	ui->setupUi(this);

	const auto setUpGroupBox = [](const auto& range, auto* minLineEdit, auto* maxLineEdit, auto* doubleSpinBox) {
		minLineEdit->setText(QString::number(range.min));
		maxLineEdit->setText(QString::number(range.max));
		// doubleSpinBox->setMinimum(min);
		// doubleSpinBox->setMaximum(max);
		doubleSpinBox->setValue((range.max + range.min) / 2.0);
	};

	setUpGroupBox(FFE_PARAMS.lambda1, ui->lambda1MinLineEdit, ui->lambda1MaxLineEdit, ui->lambda1DoubleSpinBox);
	setUpGroupBox(FFE_PARAMS.lambda2, ui->lambda2MinLineEdit, ui->lambda2MaxLineEdit, ui->lambda2DoubleSpinBox);
	setUpGroupBox(FFE_PARAMS.mu, ui->muMinLineEdit, ui->muMaxLineEdit, ui->muDoubleSpinBox);
	setUpGroupBox(FFE_PARAMS.sigma_lambda1, ui->sigmaLambda1MinLineEdit, ui->sigmaLambda1MaxLineEdit, ui->sigmaLambda1DoubleSpinBox);
	setUpGroupBox(FFE_PARAMS.sigma_lambda2, ui->sigmaLambda2MinLineEdit, ui->sigmaLambda2MaxLineEdit, ui->sigmaLambda2DoubleSpinBox);

	const auto setUpTableWidget = [](auto* tableWidget) {
		tableWidget->verticalHeader()->setVisible(false);
		tableWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
	};

	setUpTableWidget(ui->fullFactorialExperimentTableWidget);
	setUpTableWidget(ui->fractionalFactorialExperimentTableWidget);

	ui->fullFactorialExperimentTableWidget->verticalHeader()->setVisible(false);
	ui->fullFactorialExperimentTableWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);

	result = FullFactorialExperiment(FFE_PARAMS);

	const auto insertRow = [&](const FfeTableRow& row) {
		const auto rows = ui->fullFactorialExperimentTableWidget->rowCount();
		ui->fullFactorialExperimentTableWidget->insertRow(rows);
		ui->fullFactorialExperimentTableWidget->setItem(rows, 0, new QTableWidgetItem(QString::number(row.index)));
		ui->fullFactorialExperimentTableWidget->setItem(rows, 1, new QTableWidgetItem(QString::number(row.x1)));
		ui->fullFactorialExperimentTableWidget->setItem(rows, 2, new QTableWidgetItem(QString::number(row.x2)));
		ui->fullFactorialExperimentTableWidget->setItem(rows, 3, new QTableWidgetItem(QString::number(row.x3)));
		ui->fullFactorialExperimentTableWidget->setItem(rows, 4, new QTableWidgetItem(QString::number(row.x4)));
		ui->fullFactorialExperimentTableWidget->setItem(rows, 5, new QTableWidgetItem(QString::number(row.x5)));
		ui->fullFactorialExperimentTableWidget->setItem(rows, 6, new QTableWidgetItem(QString::number(row.y_mean)));
		ui->fullFactorialExperimentTableWidget->setItem(rows, 7, new QTableWidgetItem(QString::number(row.y_var)));
		ui->fullFactorialExperimentTableWidget->setItem(rows, 8, new QTableWidgetItem(QString::number(row.y_hat)));
		ui->fullFactorialExperimentTableWidget->setItem(rows, 9, new QTableWidgetItem(QString::number(row.dy_hat)));
		ui->fullFactorialExperimentTableWidget->setItem(rows, 10, new QTableWidgetItem(QString::number(row.u_hat)));
		ui->fullFactorialExperimentTableWidget->setItem(rows, 11, new QTableWidgetItem(QString::number(row.du_hat)));
	};

	for (auto&& row : result.table) {
		insertRow(row);
	}
}

MainWindow::~MainWindow()
{
	delete ui;
}

void MainWindow::on_calculatePushButton_clicked()
{
	const DotParameters dot_params{
		.lambda1       = ui->lambda1DoubleSpinBox->value(),
		.lambda2       = ui->lambda2DoubleSpinBox->value(),
		.mu            = ui->muDoubleSpinBox->value(),
		.sigma_lambda1 = ui->sigmaLambda1DoubleSpinBox->value(),
		.sigma_lambda2 = ui->sigmaLambda2DoubleSpinBox->value(),
	};

	const auto [est, act] = CalculateDot(FFE_PARAMS, result.coefficients, dot_params);

	ui->estimatedAverageWaitingTimeLineEdit->setText(QString::number(est));
	ui->actualAverageWaitingTimeLineEdit->setText(QString::number(act));
}
