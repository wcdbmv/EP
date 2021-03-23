#include "mainwindow.hpp"
#include "./ui_mainwindow.h"

#include "random.hpp"


MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
	, ui(new Ui::MainWindow)
{
	ui->setupUi(this);

	const auto setUpGroupBox = [](auto min, auto max, auto minLineEdit, auto maxLineEdit, auto doubleSpinBox) {
		minLineEdit->setText(QString::number(min));
		maxLineEdit->setText(QString::number(max));
		// doubleSpinBox->setMinimum(min);
		// doubleSpinBox->setMaximum(max);
		doubleSpinBox->setValue((max + min) / 2.0);
	};

	setUpGroupBox(LAMBDA_MIN, LAMBDA_MAX, ui->lambdaMinLineEdit, ui->lambdaMaxLineEdit, ui->lambdaDoubleSpinBox);
	setUpGroupBox(SIGMA_LAMBDA_MIN, SIGMA_LAMBDA_MAX, ui->sigmaLambdaMinLineEdit, ui->sigmaLambdaMaxLineEdit, ui->sigmaLambdaDoubleSpinBox);
	setUpGroupBox(MU_MIN, MU_MAX, ui->muMinLineEdit, ui->muMaxLineEdit, ui->muDoubleSpinBox);

	ui->ffeTableWidget->verticalHeader()->setVisible(false);
	ui->ffeTableWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);

	result = FullFactorialExperiment(PARAMS);

	const auto insertRow = [&](const FfeTableRow& row) {
		const auto rows = ui->ffeTableWidget->rowCount();
		ui->ffeTableWidget->insertRow(rows);
		ui->ffeTableWidget->setItem(rows, 0, new QTableWidgetItem(QString::number(row.index)));
		ui->ffeTableWidget->setItem(rows, 1, new QTableWidgetItem(QString::number(row.x1)));
		ui->ffeTableWidget->setItem(rows, 2, new QTableWidgetItem(QString::number(row.x2)));
		ui->ffeTableWidget->setItem(rows, 3, new QTableWidgetItem(QString::number(row.x3)));
		ui->ffeTableWidget->setItem(rows, 4, new QTableWidgetItem(QString::number(row.y_mean)));
		ui->ffeTableWidget->setItem(rows, 5, new QTableWidgetItem(QString::number(row.y_var)));
		ui->ffeTableWidget->setItem(rows, 6, new QTableWidgetItem(QString::number(row.partial_nonlinear)));
		ui->ffeTableWidget->setItem(rows, 7, new QTableWidgetItem(QString::number(row.dpn)));
		ui->ffeTableWidget->setItem(rows, 8, new QTableWidgetItem(QString::number(row.linear)));
		ui->ffeTableWidget->setItem(rows, 9, new QTableWidgetItem(QString::number(row.dl)));
	};

	for (auto&& row : result.rows) {
		insertRow(row);
	}

	ui->a0LineEdit->setText(QString::number(result.planning_matrix.a0));
	ui->a1LineEdit->setText(QString::number(result.planning_matrix.a1));
	ui->a2LineEdit->setText(QString::number(result.planning_matrix.a2));
	ui->a3LineEdit->setText(QString::number(result.planning_matrix.a3));
	ui->a12LineEdit->setText(QString::number(result.planning_matrix.a12));
	ui->a13LineEdit->setText(QString::number(result.planning_matrix.a13));
	ui->a23LineEdit->setText(QString::number(result.planning_matrix.a23));
	ui->a123LineEdit->setText(QString::number(result.planning_matrix.a123));

	ui->cochranTestLineEdit->setText(QString::number(result.cochran_test));
	ui->reproducibilityVarLineEdit->setText(QString::number(result.reproducibility_var));
	ui->adequacyVarLineEdit->setText(QString::number(result.adequacy_var));
	ui->fTestLineEdit->setText(QString::number(result.f_test));
}

MainWindow::~MainWindow()
{
	delete ui;
}

void MainWindow::on_calculatePushButton_clicked()
{
	const auto lambda = ui->lambdaDoubleSpinBox->value();
	const auto sigma_lambda = ui->sigmaLambdaDoubleSpinBox->value();
	const auto mu = ui->muDoubleSpinBox->value();

	const auto [est, act] = CalculateDot(PARAMS, result.planning_matrix, lambda, sigma_lambda, mu);

	ui->estimatedAverageWaitingTimeLineEdit->setText(QString::number(est));
	ui->actualAverageWaitingTimeLineEdit->setText(QString::number(act));
}
