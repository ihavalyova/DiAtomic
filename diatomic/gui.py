import sys
import numpy as np
from PyQt5 import QtCore, QtGui, QtWidgets
import pyqtgraph as pg

try:
    from ui_diatom import Ui_MainWindow
except (ImportError, KeyError) as exc:
    pass
    # qtCreatorFile = "%s/UI/es.ui"%lib_path
    # Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

from diatomic import *

state_count = 0
opr_count = 0

class Ui_Plot():

    def __init__(self):
        self.colors = [
            '#0066ff', '#ff0000', '#146C14', '#398989',
            '#8C36A8', '#CE9237', '#BABA55', '#1E5993',
            '#B14D6E', '#339933', '#666699'
        ]
        # self.colorDialog = QtWidgets.QColorDialog()
        # self.colorDialog.setOption(QtWidgets.QColorDialog.ShowAlphaChannel, True)
        # self.colorDialog.setOption(QtWidgets.QColorDialog.DontUseNativeDialog, True)
        # colorz = self.colorDialog.getColor(options=QtWidgets.QColorDialog.DontUseNativeDialog).name()

    def plot_eigenvectors(self, graphWidget1, rgrid, evecs_save):

        self.graphWidget1 = graphWidget1

        cind = 0
        for nlevel in range(0, evecs_save.shape[1]):
            evec = evecs_save[:, nlevel]  # substrct 1 since the counting starts from 0
            self.graphWidget1.plot(rgrid, evec, pen={'color': self.colors[cind], 'width': 1.5}, symbol=None)

            if cind > len(self.colors):
                cind = 0
            cind += 1

        styles = {'color': 'black', 'font-size': '12px'}
        self.graphWidget1.setLabel('left', 'Wavefunction', **styles)
        self.graphWidget1.setLabel('bottom', 'R (Angstrom)', **styles)
        self.graphWidget1.addLegend()

    def plot_all_radial_functions(self, graphWidget3, rgrid, ygrid):
        pass

    def plot_single_radial_functions(self, graphWidget2, xpoints, ypoints):

        self.graphWidget2 = graphWidget2

        self.curve = self.graphWidget2.plot(xpoints, ypoints, pen={'color': self.colors[0], 'width': 1.5}, symbol='o', symbolSize=7)
        styles = {'color': 'black', 'font-size': '12px'}
        self.graphWidget2.setLabel('left', 'Radial function', **styles)
        self.graphWidget2.setLabel('bottom', 'R (Angstrom)', **styles)
        self.graphWidget2.addLegend()
        self.statusBar = QtWidgets.QStatusBar()
        # self.setStatusBar(self.statusBar)
        self.curve.scene().sigMouseMoved.connect(self.onMouseMoved)

    def onMouseMoved(self, point):
        p = self.graphWidget2.plotItem.vb.mapSceneToView(point)
        self.statusBar.showMessage("{}-{}".format(p.x(), p.y()))

class Diatom(QtWidgets.QMainWindow, Ui_MainWindow):

    def __init__(self):

        super().__init__()
        Ui_MainWindow.__init__(self)
        self.setupUi(self)

        self.ui_plot = Ui_Plot()
        self.graphWidget1.setBackground('w')
        self.graphWidget2.setBackground('w')

        self.sel_isotopes = None
        self.ref_energy = None
        self.refj = None
        self.input_params_file = None
        self.obs_enr_file = None
        self.ngrid = 0
        self.rmin = 0.0
        self.rmax = 0.0
        self.solver = 'sinc'
        self.state_labels = [1]
        self.jstart = [0.0]
        self.jend = [0.0]
        self.jvalues = [[]]
        self.lambdas = [0]
        self.spins = [0.0]
        self.sigmas = [0.0]
        self.esymm = [True]
        self.fsymm = [False]
        self.save_eigenvec = False
        self.nlevels_save = [0]

        self.opr_type_names = [
            'Potential Energy',
            'Spin-orbit',
            'L-uncoupling',
            'Spin-uncoupling',
            'Spin-rotation',
            'Spin-spin'
        ]
        self.model_types = [
            'Pointwise',
            'Morse',
        ]

        self.opr_labels = [1]
        self.opr_types = [self.opr_type_names[0]]
        self.pair_state1 = [1]
        self.pair_state2 = [1]
        self.opr_models = [self.model_types[0]]
        self.opr_rad_labels = ['']
        self.opr_rot_corr = [0.0]
        self.enr_range1 = None
        self.enr_range2 = None
        self.enr_index1 = None
        self.enr_index2 = None
        self.calc_eigenenr = True
        self.calc_enr_levels = False
        self.ident_opt = 1

        self._init_state_dynamic_widget_lists()
        self._init_operator_dynamic_widget_lists()

        self.moleculeLineEdit.editingFinished.connect(self._get_molecule)
        self.refJLineEdit.editingFinished.connect(self._get_refj)
        self.refEnergyLineEdit.editingFinished.connect(self._get_ref_energy)

        self._toggle_table_widget()

        self.paramsFileToolButton.clicked.connect(
            self._on_input_parameters_file_button_clicked
        )
        self.previewParamsPushButton.clicked.connect(
            self._on_preview_input_params_button_clicked
        )
        self.obsEnergyFileToolButton.clicked.connect(
            self._on_obs_energies_file_button_clicked
        )
        self.previewObsEnrPushButton.clicked.connect(
            self._on_preview_obs_enr_button_clicked
        )

        self.nGridPointsLineEdit.editingFinished.connect(self._get_ngrid_points)
        self.rminLineEdit.editingFinished.connect(self._get_rmin_grid)
        self.rmaxLineEdit.editingFinished.connect(self._get_rmax_grid)

        self.noneqidGridRadioButton.clicked.connect(
            self._on_noneqid_grid_radio_button_clicked
        )
        self.sincSolverRadioButton.clicked.connect(
            self._on_sinc_solver_radio_button_clicked
        )
        self.fourierSolverRadioButton.clicked.connect(
            self._on_fourier_solver_radio_button_clicked
        )
        self.fd5SolverRadioButton.clicked.connect(
            self._on_fd5_solver_radio_button_clicked
        )

        self.addStatePushButton.clicked.connect(
            self._on_add_state_button_clicked
        )

        self.stateJsSpinBox[0].valueChanged.connect(
            lambda: self._on_state_js_spin_box_value_change(0)
        )
        self.stateJeSpinBox[0].valueChanged.connect(
            lambda: self._on_state_je_spin_box_value_change(0)
        )
        self.stateJvLineEdit[0].editingFinished.connect(
            lambda: self._get_jvalues(0)
        )
        self.stateLambdaSpinBox[0].valueChanged.connect(
            lambda: self._on_state_lambda_spin_box_value_change(0)
        )
        self.stateSpinSpinBox[0].valueChanged.connect(
            lambda: self._on_state_spin_spin_box_value_change(0)
        )
        self.stateSigmaSpinBox[0].valueChanged.connect(
            lambda: self._on_state_sigma_spin_box_value_change(0)
        )
        self.stateECheckBox[0].stateChanged.connect(
            lambda: self._on_state_esymm_state_change(0)
        )
        self.stateFCheckBox[0].stateChanged.connect(
            lambda: self._on_state_fsymm_state_change(0)
        )

        # self.removeStatePushButton.clicked.connect(
        #     self._on_remove_state_button_clicked
        # )

        self.operatorTypeComboBox1.addItems(self.opr_type_names)
        self.modelTypeComboBox1.addItems(self.model_types)

        self.addOperatorPushButton.clicked.connect(
            self._on_add_operator_button_clicked
        )

        self.operatorTypeComboBox[0].currentTextChanged.connect(
            lambda: self._on_operator_type_change(0)
        )
        self.modelTypeComboBox[0].currentTextChanged.connect(
            lambda: self._on_operator_model_change(0)
        )
        self.pairStates1SpinBox[0].valueChanged.connect(
            lambda: self._on_operator_pair_state1_value_change(0)
        )
        self.pairStates2SpinBox[0].valueChanged.connect(
            lambda: self._on_operator_pair_state2_value_change(0)
        )
        self.labelRadParamsLineEdit[0].editingFinished.connect(
            lambda: self._get_operator_label(0)
        )
        self.labelRotCorrLineEdit[0].editingFinished.connect(
            lambda: self._get_operator_rot_corr(0)
        )
        self.showRadParamsPushButton[0].clicked.connect(
           lambda: self._on_show_radial_parameters_button_clicked(0)
        )
        self.plotRadParamsPushButton1.clicked.connect(
            lambda: self._on_plot_rad_params_push_button_clicked(0)
        )

        self.energyRange1LineEdit.editingFinished.connect(self._get_energy_range1)
        self.energyRange2LineEdit.editingFinished.connect(self._get_energy_range2)
        self.energyIndex1LineEdit.editingFinished.connect(self._get_energy_index1)
        self.energyIndex2LineEdit.editingFinished.connect(self._get_energy_index2)

        self.calcEigenenergyCheckBox.stateChanged.connect(
            self._on_calc_eigenenr_state_change
        )
        self.calcEnergyLevelsCheckBox.stateChanged.connect(
            self._on_calc_enr_levels_state_change
        )
        self.identLevelsByContrRadioButton.clicked.connect(
            self._on_ident_levels_by_contr_button_clicked
        )
        self.identLevelsByEnrRadioButton.clicked.connect(
            self._on_ident_levels_by_enr_button_clicked
        )

        self.saveEigenvectrosCheckBox.stateChanged.connect(
            self._on_save_eigenvectors_state_change
        )
        self.plotNlevelsLineEdit.editingFinished.connect(
            self._get_nlevels_to_save
        )
        self.previewOutputPushButton.clicked.connect(
            self._on_preview_output_button_clicked
        )

        self.solvePushButton.clicked.connect(
            self._on_solve_button_clicked
        )

        self.svdOptMethodRadioButton.clicked.connect(
            self._on_svd_opt_method_button_clicked
        )
        self.minuitOptMethodRadioButton.clicked.connect(
            self._on_minuit_opt_method_button_clicked
        )
        self.toleranceLineEdit.editingFinished.connect(
            self._get_tolerance_value
        )
        self.numericalDerivRadioButton.clicked.connect(
            self._on_numerical_deriv_button_clicked
        )
        self.analyticalDerivRadioButton.clicked.connect(
            self._on_analytical_deriv_button_clicked
        )
        self.numberOfIterLineEdit.editingFinished.connect(
            self._get_number_of_iterations
        )

    def _get_molecule(self):

        self.molecule = self.moleculeLineEdit.text()
        print(self.molecule)

        # list all possible isotopes and get the selected ones
        self.isotopes, self.atoms_info = Diatomic.get_atomic_isotopes(self.molecule)
        print(self.isotopes)

        self.isotopesListWidget.clear()

        self.isotopesListWidget.setSelectionMode(
            QtWidgets.QAbstractItemView.ExtendedSelection
        )
        for _, isotope in enumerate(self.isotopes):
            item = QtWidgets.QListWidgetItem(isotope)
            self.isotopesListWidget.addItem(item)
        self.isotopesListWidget.itemClicked.connect(self._get_selected_isotopes)

    def _get_selected_isotopes(self):

        items = self.isotopesListWidget.selectedItems()
        self.sel_isotopes = []
        sel_iso_inds = []
        for i in range(len(items)):
            sel_iso = str(self.isotopesListWidget.selectedItems()[i].text())
            self.sel_isotopes.append(sel_iso)
            sel_iso_inds.append(self.isotopes.index(sel_iso))

        sel_iso_text = 'Selected isotopologues: '
        for i in range(len(self.sel_isotopes)):
            sel_iso_text += f'{self.sel_isotopes[i]} '

        sel_iso_info = []
        for sel_iso_ind in sel_iso_inds:
            sel_iso_info.append(self.atoms_info[sel_iso_ind])

        self.isotopesTextEdit.setPlainText(sel_iso_text + '\n\n' + ''.join(sel_iso_info))

    def _get_refj(self):
        self.refj = self.refJLineEdit.text()
        print(self.refj)

    def _get_ref_energy(self):
        self.ref_energy = self.refEnergyLineEdit.text()
        print(self.ref_energy)

    def _toggle_table_widget(self):
        self.tableWidget.setHidden(not self.tableWidget.isHidden())

    def _on_input_parameters_file_button_clicked(self):

        self.input_params_file, _ = QtWidgets.QFileDialog.getOpenFileName(
            parent=self, caption='Open file', directory='.', filter='*')
        if self.input_params_file:
            self.paramsFileLineEdit.setText(self.input_params_file)

    def _on_preview_input_params_button_clicked(self):

        self._toggle_table_widget()

        params = []
        with open(self.input_params_file, encoding='utf-8') as fstrm:
            for line in fstrm:
                if not line.startswith('@'):
                    params.append(line.split())

        nrows = len(params)
        ncols = max([len(p) for p in params])
        self.tableWidget.setRowCount(nrows)
        self.tableWidget.setColumnCount(ncols)

        for irow, rvalue in enumerate(params):
            for icol, cvalue in enumerate(rvalue):
                new_item = QtWidgets.QTableWidgetItem(cvalue)
                self.tableWidget.setItem(irow, icol, new_item)

    def _on_obs_energies_file_button_clicked(self):

        self.obs_enr_file, _ = QtWidgets.QFileDialog.getOpenFileName(
            parent=self, caption='Open file', directory='.', filter='*')
        if self.obs_enr_file:
            self.obsEnergyFileLineEdit.setText(self.obs_enr_file)

    def _on_preview_obs_enr_button_clicked(self):

        try:
            obs_data = np.genfromtxt(self.obs_enr_file, comments='#', autostrip=True)
        except Exception as excp:
            raise SystemExit(excp) from excp

        nrows, ncols = obs_data.shape
        self.tableWidget_2.setRowCount(nrows)
        self.tableWidget_2.setColumnCount(ncols)

        for irow, rvalue in enumerate(obs_data):
            for icol, cvalue in enumerate(rvalue):
                cvalue_str = '{:10.5f}'.format(cvalue)
                new_item = QtWidgets.QTableWidgetItem(cvalue_str)
                self.tableWidget_2.setItem(irow, icol, new_item)

    def _on_preview_output_button_clicked(self):

        if self.calc_eigenenr:
            try:
                cal_eignenr = np.genfromtxt(self.file_eigenenr, comments='#', autostrip=True)
            except Exception as excp:
                raise SystemExit(excp) from excp

            nrows, ncols = cal_eignenr.shape
            self.tableWidget_5.setRowCount(nrows)
            self.tableWidget_5.setColumnCount(ncols)

            for irow, rvalue in enumerate(cal_eignenr):
                for icol, cvalue in enumerate(rvalue):
                    cvalue_str = '{:10.5f}'.format(cvalue)
                    new_item = QtWidgets.QTableWidgetItem(cvalue_str)
                    self.tableWidget_5.setItem(irow, icol, new_item)

        if self.calc_enr_levels:
            try:
                cal_enr_levs = np.genfromtxt(self.file_enr_levels, comments='#', autostrip=True)
            except Exception as excp:
                raise SystemExit(excp) from excp

            nrows, ncols = cal_enr_levs.shape
            self.tableWidget_6.setRowCount(nrows)
            self.tableWidget_6.setColumnCount(ncols)

            for irow, rvalue in enumerate(cal_enr_levs):
                for icol, cvalue in enumerate(rvalue):
                    cvalue_str = '{:10.5f}'.format(cvalue)
                    new_item = QtWidgets.QTableWidgetItem(cvalue_str)
                    self.tableWidget_6.setItem(irow, icol, new_item)

        if self.save_eigenvec:
            try:
                eigenvecs = np.genfromtxt(self.file_eigenvecs, comments='#', autostrip=True)
            except Exception as excp:
                raise SystemExit(excp) from excp

            nrows, ncols = eigenvecs.shape
            self.tableWidget_7.setRowCount(nrows)
            self.tableWidget_7.setColumnCount(ncols)

            for irow, rvalue in enumerate(eigenvecs):
                for icol, cvalue in enumerate(rvalue):
                    cvalue_str = '{:10.5f}'.format(cvalue)
                    new_item = QtWidgets.QTableWidgetItem(cvalue_str)
                    self.tableWidget_7.setItem(irow, icol, new_item)

    def _get_ngrid_points(self):
        self.ngrid = int(self.nGridPointsLineEdit.text())

    def _get_rmin_grid(self):
        self.rmin = float(self.rminLineEdit.text())

    def _get_rmax_grid(self):
        self.rmax = float(self.rmaxLineEdit.text())

    def _on_noneqid_grid_radio_button_clicked(self):

        if self.noneqidGridRadioButton.isChecked():
            msg = QtWidgets.QMessageBox()
            msg.setIcon(QtWidgets.QMessageBox.Information)
            msg.setText("The non-equdistant calculations are currently not supported!")
            # msg.setInformativeText('More information')
            msg.setWindowTitle("Information")
            msg.exec_()

    def _on_sinc_solver_radio_button_clicked(self):
        self.solver = 'sinc'

    def _on_fd5_solver_radio_button_clicked(self):
        self.solver = 'fourier'

    def _on_fourier_solver_radio_button_clicked(self):
        self.solver = 'fd5'

    def _on_add_state_button_clicked(self):
        self._generate_state_widgets()

    def _on_add_operator_button_clicked(self):
        self._generate_operator_widgets()

    def _on_remove_state_button_clicked(self):
        self._remove_state_widgets()

    def _on_state_js_spin_box_value_change(self, index):
        self.jstart[index] = self.stateJsSpinBox[index].value()

    def _on_state_je_spin_box_value_change(self, index):
        self.jend[index] = self.stateJeSpinBox[index].value()

    def _get_jvalues(self, index):
        jvalues_str = self.stateJvLineEdit[index].text().split()
        self.jvalues[index] = [float(i) for i in jvalues_str]

    def _on_state_esymm_state_change(self, index):

        if self.stateECheckBox[index].isChecked():
            self.esymm[index] = True
        else:
            self.esymm[index] = False

    def _on_state_fsymm_state_change(self, index):

        if self.stateFCheckBox[index].isChecked():
            self.fsymm[index] = True
        else:
            self.fsymm[index] = False

    def _on_state_lambda_spin_box_value_change(self, index):
        self.lambdas[index] = self.stateLambdaSpinBox[index].value()

    def _on_state_spin_spin_box_value_change(self, index):
        self.spins[index] = self.stateSpinSpinBox[index].value()

    def _on_state_sigma_spin_box_value_change(self, index):
        self.sigmas[index] = self.stateSigmaSpinBox[index].value()

    def _on_operator_type_change(self, index):
        self.opr_types[index] = self.operatorTypeComboBox[index].currentText()

    def _on_operator_model_change(self, index):
        self.opr_models[index] = self.modelTypeComboBox[index].currentText()

    def _on_operator_pair_state1_value_change(self, index):
        self.pair_state1[index] = self.pairStates1SpinBox[index].value()

    def _on_operator_pair_state2_value_change(self, index):
        self.pair_state2[index] = self.pairStates2SpinBox[index].value()

    def _get_operator_label(self, index):
        self.opr_rad_labels[index] = self.labelRadParamsLineEdit[index].text()

    def _get_operator_rot_corr(self, index):
        self.opr_rot_corr[index] = float(self.labelRotCorrLineEdit[index].text())

    def _on_show_radial_parameters_button_clicked(self, index):

        # self._toggle_table_widget()

        params = []
        found = False
        with open(self.input_params_file, encoding='utf-8') as fstrm:
            for line in fstrm:
                if line.startswith('@'+self.opr_rad_labels[index]+':'):
                    found = True
                    continue
                elif line.startswith('@') and found:
                    break

                if found:
                    params.append(line.split())

        nrows = len(params)
        ncols = max([len(p) for p in params])
        self.tableWidget.setRowCount(nrows)
        self.tableWidget.setColumnCount(ncols)

        for irow, rvalue in enumerate(params):
            for icol, cvalue in enumerate(rvalue):
                new_item = QtWidgets.QTableWidgetItem(cvalue)
                self.tableWidget.setItem(irow, icol, new_item)

    def _get_energy_range1(self):
        self.enr_range1 = float(self.energyRange1LineEdit.text())

    def _get_energy_range2(self):
        self.enr_range2 = float(self.energyRange2LineEdit.text())

    def _get_energy_index1(self):
        self.enr_index1 = int(self.energyIndex1LineEdit.text())

    def _get_energy_index2(self):
        self.enr_index2 = int(self.energyIndex2LineEdit.text())

    def _on_calc_eigenenr_state_change(self):

        if self.calcEigenenergyCheckBox.isChecked():
            self.calc_eigenenr = True
        else:
            self.calc_eigenenr = False

    def _on_calc_enr_levels_state_change(self):

        if self.calcEnergyLevelsCheckBox.isChecked():
            self.calc_enr_levels = True
        else:
            self.calc_enr_levels = False

    def _on_ident_levels_by_contr_button_clicked(self):
        self.ident_opt = 1

    def _on_ident_levels_by_enr_button_clicked(self):
        self.ident_opt = 0

    def _get_ef_symmetry(self, istate):

        if self.esymm[istate] and self.fsymm[istate]:
            return (0, 1)
        if self.esymm[istate] and not self.fsymm[istate]:
            return (1,)

        return (0,)

    def _on_save_eigenvectors_state_change(self):

        if self.saveEigenvectrosCheckBox.isChecked():
            self.save_eigenvec = True
        else:
            self.save_eigenvec = False

    def _get_nlevels_to_save(self):
        nlevels_save_str = self.plotNlevelsLineEdit.text().split()
        self.nlevels_save = [int(i) for i in nlevels_save_str]

    def _on_plot_rad_params_push_button_clicked(self, index):

        params = []
        found = False
        with open(self.input_params_file, encoding='utf-8') as fstrm:
            for line in fstrm:
                if line.startswith('@'+self.opr_rad_labels[index]+':'):
                    found = True
                    continue
                elif line.startswith('@') and found:
                    break

                if found:
                    params.append(line.split())

        params_arr = np.array(params, dtype=np.float64)

        self.ui_plot.plot_single_radial_functions(self.graphWidget2, params_arr[:, 0], params_arr[:, 1])

    def _on_svd_opt_method_button_clicked(self):
        pass

    def _on_minuit_opt_method_button_clicked(self):
        pass

    def _get_tolerance_value(self):
        pass

    def _get_number_of_iterations(self):
        pass

    def _on_numerical_deriv_button_clicked(self):
        pass

    def _on_analytical_deriv_button_clicked(self):
        pass

    def _on_solve_button_clicked(self):

        diatomic = Diatomic(
            molecule=self.sel_isotopes,
            refj=self.refj,
            ref_enr=self.ref_energy
        )

        diatomic.set_data_parameters(self.input_params_file)

        if self.obs_enr_file is not None:
            diatomic.set_observed_energies(self.obs_enr_file)

        grid = Grid(
            npoints=self.ngrid,
            rgrid=(self.rmin, self.rmax),
            solver=self.solver
        )

        states = []
        for _, istate in enumerate(self.state_labels):
            states.append(Basis(
                istate,
                Js=self.jstart[istate-1],
                Je=self.jend[istate-1],
                symmetry=self._get_ef_symmetry(istate-1),
                _lambda=self.lambdas[istate-1],
                spin=self.spins[istate-1],
                sigma=self.sigmas[istate-1]
            ))

        objs = (diatomic, grid, states)

        operators = []
        for _, iopr in enumerate(self.opr_labels):

            if self.opr_types[iopr-1] == self.opr_type_names[0]:
                operators.append(PotEnr(
                    objs,
                    pair_states=(self.pair_state1[iopr-1], self.pair_state2[iopr-1]),
                    model=self.opr_models[iopr-1],
                    label=self.opr_rad_labels[iopr-1],
                    rotc=self.opr_rot_corr[iopr-1]
                ))
            elif self.opr_types[iopr-1] == self.opr_type_names[1]:
                operators.append(SO(
                    objs,
                    pair_states=(self.pair_state1[iopr-1], self.pair_state2[iopr-1]),
                    model=self.opr_models[iopr-1],
                    label=self.opr_rad_labels[iopr-1],
                    rotc=self.opr_rot_corr[iopr-1]
                ))
            elif self.opr_types[iopr-1] == self.opr_type_names[2]:
                operators.append(LJ(
                    objs,
                    pair_states=(self.pair_state1[iopr-1], self.pair_state2[iopr-1]),
                    model=self.opr_models[iopr-1],
                    label=self.opr_rad_labels[iopr-1],
                    rotc=self.opr_rot_corr[iopr-1]
                ))
            elif self.opr_types[iopr-1] == self.opr_type_names[3]:
                operators.append(PotEnr(
                    SJ,
                    pair_states=(self.pair_state1[iopr-1], self.pair_state2[iopr-1]),
                    model=self.opr_models[iopr-1],
                    label=self.opr_rad_labels[iopr-1],
                    rotc=self.opr_rot_corr[iopr-1]
                ))
            elif self.opr_types[iopr-1] == self.opr_type_names[4]:
                operators.append(SR(
                    objs,
                    pair_states=(self.pair_state1[iopr-1], self.pair_state2[iopr-1]),
                    model=self.opr_models[iopr-1],
                    label=self.opr_rad_labels[iopr-1],
                    rotc=self.opr_rot_corr[iopr-1]
                ))
            elif self.opr_types[iopr-1] == self.opr_type_names[5]:
                operators.append(SS(
                    objs,
                    pair_states=(self.pair_state1[iopr-1], self.pair_state2[iopr-1]),
                    model=self.opr_models[iopr-1],
                    label=self.opr_rad_labels[iopr-1],
                    rotc=self.opr_rot_corr[iopr-1]
                ))

        H = Hamiltonian(objs)

        self.rgrid = H.rgrid

        energies, _ = H.solve(
            operators,
            energy_range_value=(self.enr_range1, self.enr_range2)
        )

        if self.calc_eigenenr:
            self.file_eigenenr = f'{self.molecule}_eigenenergies.dat'
            H.save_full_energy_list(energies, filename=self.file_eigenenr)

        if self.calc_enr_levels:
            levels = H.get_energy_levels(ident_option=self.ident_opt)
            self.file_enr_levels = f'{self.molecule}_energy_levels.dat'
            H.save_energy_list(filename=self.file_enr_levels)

        if self.save_eigenvec:
            self.evecs_matrix = H.evecs_matrix

            # review the list of eigenvectors
            nlevels_save_rev = []
            for nlevel in self.nlevels_save:
                if nlevel < self.evecs_matrix.shape[1]:
                    nlevels_save_rev.append(nlevel)

            # store eigenvectors
            try:
                evecs_save = self.evecs_matrix[:, self.nlevels_save]
            except IndexError:
                evecs_save = self.evecs_matrix[:, nlevels_save_rev]
            self.file_eigenvecs = f'{self.molecule}_eigenvectors.dat'
            np.savetxt(self.file_eigenvecs, evecs_save)

            # plot eigenvectors
            self.ui_plot.plot_eigenvectors(self.graphWidget1, self.rgrid, evecs_save)

        print(self.state_labels)
        print(self.jstart)
        print(self.jend)
        print(self.jvalues)
        print(self.lambdas)
        print(self.spins)
        print(self.sigmas)
        print(self.esymm)
        print(self.fsymm)
        print(self.opr_types)
        print(self.opr_models)
        print(self.pair_state1)
        print(self.pair_state2)
        print(self.opr_rad_labels)
        print(self.opr_rot_corr)
        print(self.enr_range1)
        print(self.enr_range2)
        print(self.enr_index1)
        print(self.enr_index2)
        print(self.calc_eigenenr)
        print(self.calc_enr_levels)
        print(self.ident_opt)

    def _init_state_dynamic_widget_lists(self):

        self.stateGridLayouts = []
        self.stateLabels = []
        self.stateJsLabels = []
        self.stateJeLabels = []
        self.stateJvLabels = []
        self.stateLabelSpinBox = []
        self.stateJsSpinBox = []
        self.stateJeSpinBox = []
        self.stateJvLineEdit = []
        self.stateSymmetryLabels = []
        self.stateLambdaLabels = []
        self.stateSpinLabels = []
        self.stateSigmaLabels = []
        self.stateECheckBox = []
        self.stateFCheckBox = []
        self.stateLambdaSpinBox = []
        self.stateSpinSpinBox = []
        self.stateSigmaSpinBox = []

        self.stateGridLayouts.append(self.gridLayout_4)
        self.stateLabels.append(self.label)
        self.stateJsLabels.append(self.label_2)
        self.stateJeLabels.append(self.label_3)
        self.stateJvLabels.append(self.label_9)
        self.stateLabelSpinBox.append(self.spinBox_4)
        self.stateJsSpinBox.append(self.doubleSpinBox_3)
        self.stateJeSpinBox.append(self.doubleSpinBox_4)
        self.stateJvLineEdit.append(self.lineEdit_2)
        self.stateSymmetryLabels.append(self.label_4)
        self.stateLambdaLabels.append(self.label_6)
        self.stateSpinLabels.append(self.label_7)
        self.stateSigmaLabels.append(self.label_8)
        self.stateECheckBox.append(self.checkBox)
        self.stateFCheckBox.append(self.checkBox_2)
        self.stateLambdaSpinBox.append(self.spinBox_3)
        self.stateSpinSpinBox.append(self.doubleSpinBox)
        self.stateSigmaSpinBox.append(self.doubleSpinBox_2)

    def _remove_state_widgets(self):
        self._remove_layout(self.stateGridLayouts)

    def _remove_layout(self, layout):
        print(layout)
        for i in range(len(layout)):
            print(i, layout[i])
            temp_layout = layout[i]
            if temp_layout is not None:
                widget = temp_layout.widget()
                if widget is not None:
                    widget.deleteLater()
            else:
                return
            if temp_layout.layout() is not None:
                self._remove_layout(temp_layout.layout())

    def _generate_state_widgets(self):

        global state_count
        state_count += 1
        top_margin = 30

        self.jstart.append(0.0)
        self.jend.append(0.0)
        self.jvalues.append([])
        self.lambdas.append(0)
        self.spins.append(0.0)
        self.sigmas.append(0.0)
        self.esymm.append(True)
        self.fsymm.append(False)

        self.stateGridLayouts.append(QtWidgets.QGridLayout())
        self.stateGridLayouts[state_count].setObjectName(f'sgridLayout_{state_count}')
        self.stateGridLayouts[state_count].setContentsMargins(0, top_margin, 0, 0)

        self.stateLabels.append(QtWidgets.QLabel(self.scrollAreaWidgetContents))
        self.stateLabels[state_count].setObjectName(f'slabel_{state_count}')
        self.stateGridLayouts[state_count].addWidget(self.stateLabels[state_count], 0, 0, 1, 1)
        self.stateLabels[state_count].setText("State Label:")
        state_label = state_count+1
        self.state_labels.append(state_label)

        self.stateJsLabels.append(QtWidgets.QLabel(self.scrollAreaWidgetContents))
        self.stateJsLabels[state_count].setObjectName(f'jslabel_{state_count}')
        self.stateGridLayouts[state_count].addWidget(self.stateJsLabels[state_count], 1, 0, 1, 1)
        self.stateJsLabels[state_count].setText("J start =")

        self.stateJeLabels.append(QtWidgets.QLabel(self.scrollAreaWidgetContents))
        self.stateJeLabels[state_count].setObjectName(f'jelabel_{state_count}')
        self.stateGridLayouts[state_count].addWidget(self.stateJeLabels[state_count], 2, 0, 1, 1)
        self.stateJeLabels[state_count].setText("J end =")

        self.stateJvLabels.append(QtWidgets.QLabel(self.scrollAreaWidgetContents))
        self.stateJvLabels[state_count].setObjectName(f'jvlabel_{state_count}')
        self.stateGridLayouts[state_count].addWidget(self.stateJvLabels[state_count], 3, 0, 1, 1)
        self.stateJvLabels[state_count].setText("J values:")

        self.stateLabelSpinBox.append(QtWidgets.QSpinBox(self.scrollAreaWidgetContents))
        self.stateLabelSpinBox[state_count].setObjectName(f'stateLabelSpinBox_{state_count}')
        self.stateLabelSpinBox[state_count].setMinimum(state_count+1)
        self.stateGridLayouts[state_count].addWidget(self.stateLabelSpinBox[state_count], 0, 1, 1, 1)
        self.stateLabelSpinBox[state_count].setReadOnly(True)

        self.stateJsSpinBox.append(QtWidgets.QDoubleSpinBox(self.scrollAreaWidgetContents))
        self.stateJsSpinBox[state_count].setObjectName(f'JsdoubleSpinBox_{state_count}')
        self.stateJsSpinBox[state_count].setDecimals(1)
        self.stateGridLayouts[state_count].addWidget(self.stateJsSpinBox[state_count], 1, 1, 1, 1)
        self.stateJsSpinBox[state_label-1].valueChanged.connect(
            lambda: self._on_state_js_spin_box_value_change(state_label-1)
        )

        self.stateJeSpinBox.append(QtWidgets.QDoubleSpinBox(self.scrollAreaWidgetContents))
        self.stateJeSpinBox[state_count].setObjectName(f'JedoubleSpinBox_{state_count}')
        self.stateJeSpinBox[state_count].setDecimals(1)
        self.stateGridLayouts[state_count].addWidget(self.stateJeSpinBox[state_count], 2, 1, 1, 1)
        self.stateJeSpinBox[state_label-1].valueChanged.connect(
            lambda: self._on_state_je_spin_box_value_change(state_label-1)
        )

        self.stateJvLineEdit.append(QtWidgets.QLineEdit(self.scrollAreaWidgetContents))
        self.stateJvLineEdit[state_count].setObjectName(f'stateJvlineEdit_{state_count}')
        self.stateGridLayouts[state_count].addWidget(self.stateJvLineEdit[state_count], 3, 1, 1, 2)
        self.stateJvLineEdit[state_label-1].editingFinished.connect(
            lambda: self._get_jvalues(state_label-1)
        )

        self.stateSymmetryLabels.append(QtWidgets.QLabel(self.scrollAreaWidgetContents))
        self.stateSymmetryLabels[state_count].setObjectName(f'stateSymmetryLabels_{state_count}')
        self.stateGridLayouts[state_count].addWidget(self.stateSymmetryLabels[state_count], 0, 3, 1, 2)
        self.stateSymmetryLabels[state_count].setText("Symmetry:")

        self.stateLambdaLabels.append(QtWidgets.QLabel(self.scrollAreaWidgetContents))
        self.stateLambdaLabels[state_count].setObjectName(f'stateLambdaLabels_{state_count}')
        self.stateGridLayouts[state_count].addWidget(self.stateLambdaLabels[state_count], 1, 3, 1, 2)
        self.stateLambdaLabels[state_count].setText("Lambda =")

        self.stateSpinLabels.append(QtWidgets.QLabel(self.scrollAreaWidgetContents))
        self.stateSpinLabels[state_count].setObjectName(f'stateSpinLabels{state_count}')
        self.stateGridLayouts[state_count].addWidget(self.stateSpinLabels[state_count], 2, 3, 1, 2)
        self.stateSpinLabels[state_count].setText("Spin =")

        self.stateSigmaLabels.append(QtWidgets.QLabel(self.scrollAreaWidgetContents))
        self.stateSigmaLabels[state_count].setObjectName(f'stateSigmaLabels{state_count}')
        self.stateGridLayouts[state_count].addWidget(self.stateSigmaLabels[state_count], 3, 3, 1, 2)
        self.stateSigmaLabels[state_count].setText("Sigma =")

        self.stateECheckBox.append(QtWidgets.QCheckBox(self.scrollAreaWidgetContents))
        self.stateECheckBox[state_count].setObjectName("stateECheckBox")
        self.stateECheckBox[state_count].setChecked(True)
        self.stateGridLayouts[state_count].addWidget(self.stateECheckBox[state_count], 0, 5, 1, 1)
        self.stateECheckBox[state_count].setText("e")
        self.stateECheckBox[state_label-1].stateChanged.connect(
            lambda: self._on_state_esymm_state_change(state_label-1)
        )

        self.stateFCheckBox.append(QtWidgets.QCheckBox(self.scrollAreaWidgetContents))
        self.stateFCheckBox[state_count].setObjectName("stateFCheckBox")
        self.stateGridLayouts[state_count].addWidget(self.stateFCheckBox[state_count], 0, 6, 1, 1)
        self.stateFCheckBox[state_count].setText("f")
        self.stateFCheckBox[state_label-1].stateChanged.connect(
            lambda: self._on_state_fsymm_state_change(state_label-1)
        )

        self.stateLambdaSpinBox.append(QtWidgets.QSpinBox(self.scrollAreaWidgetContents))
        self.stateLambdaSpinBox[state_count].setObjectName(f'stateLambdaSpinBox{state_count}')
        self.stateLambdaSpinBox[state_count].setMinimum(0)
        self.stateGridLayouts[state_count].addWidget(self.stateLambdaSpinBox[state_count], 1, 5, 1, 1)
        self.stateLambdaSpinBox[state_label-1].valueChanged.connect(
            lambda: self._on_state_lambda_spin_box_value_change(state_label-1)
        )

        self.stateSpinSpinBox.append(QtWidgets.QDoubleSpinBox(self.scrollAreaWidgetContents))
        self.stateSpinSpinBox[state_count].setObjectName(f'stateSpinSpinBox{state_count}')
        self.stateSpinSpinBox[state_count].setDecimals(1)
        self.stateGridLayouts[state_count].addWidget(self.stateSpinSpinBox[state_count], 2, 5, 1, 1)
        self.stateSpinSpinBox[state_label-1].valueChanged.connect(
            lambda: self._on_state_spin_spin_box_value_change(state_label-1)
        )

        self.stateSigmaSpinBox.append(QtWidgets.QDoubleSpinBox(self.scrollAreaWidgetContents))
        self.stateSigmaSpinBox[state_count].setObjectName(f'stateSigmaSpinBox{state_count}')
        self.stateSigmaSpinBox[state_count].setDecimals(1)
        self.stateGridLayouts[state_count].addWidget(self.stateSigmaSpinBox[state_count], 3, 5, 1, 1)
        self.stateSigmaSpinBox[state_label-1].valueChanged.connect(
            lambda: self._on_state_sigma_spin_box_value_change(state_label-1)
        )

        self.gridLayout_7.addLayout(self.stateGridLayouts[state_count], state_count+1, 0, 1, 1)
        self.stateLabels[state_count].show()

    def _init_operator_dynamic_widget_lists(self):

        self.operatorGridLayouts = []
        self.operatorTypeLabels = []
        self.operatorPairStateLabels = []
        self.operatorModelTypeLabels = []
        self.operatorRadialLabelLabels = []
        self.operatorRotCorrLabels = []
        self.operatorTypeComboBox = []
        self.pairStates1SpinBox = []
        self.pairStates2SpinBox = []
        self.modelTypeComboBox = []
        self.labelRadParamsLineEdit = []
        self.labelRotCorrLineEdit = []
        self.showRadParamsPushButton = []
        self.plotRadParamsPushButton = []

        self.operatorGridLayouts.append(self.gridLayout_6)
        self.operatorTypeLabels.append(self.operatorTypeLabel1)
        self.operatorPairStateLabels.append(self.pairStatesLabel1)
        self.operatorModelTypeLabels.append(self.modelTypeLabel1)
        self.operatorRadialLabelLabels.append(self.labelParamsLabel1)
        self.operatorRotCorrLabels.append(self.rotCorrLabel1)
        self.operatorTypeComboBox.append(self.operatorTypeComboBox1)
        self.pairStates1SpinBox.append(self.pairStates1SpinBox1)
        self.pairStates2SpinBox.append(self.pairStates2SpinBox1)
        self.modelTypeComboBox.append(self.modelTypeComboBox1)
        self.labelRadParamsLineEdit.append(self.labelRadParamsLineEdit1)
        self.labelRotCorrLineEdit.append(self.labelRotCorrLineEdit1)
        self.showRadParamsPushButton.append(self.showRadParamsPushButton1)
        self.plotRadParamsPushButton.append(self.plotRadParamsPushButton1)

    def _generate_operator_widgets(self):

        global opr_count
        opr_count += 1
        top_margin = 30
        opr_label = opr_count+1

        self.opr_types.append(self.opr_type_names[0])
        self.opr_models.append(self.model_types[0])
        self.pair_state1.append(1)
        self.pair_state2.append(1)
        self.opr_rad_labels.append('')
        self.opr_rot_corr.append(0.0)

        self.opr_labels.append(opr_label+1)

        self.operatorGridLayouts.append(QtWidgets.QGridLayout())
        self.operatorGridLayouts[opr_count].setObjectName(f'ogridLayout_{opr_count}')
        self.operatorGridLayouts[opr_count].setContentsMargins(0, top_margin, 0, 0)

        self.operatorTypeLabels.append(QtWidgets.QLabel(self.scrollAreaWidgetContents_2))
        self.operatorTypeLabels[opr_count].setObjectName(f'ooperatorTypeLabel_{opr_count}')
        self.operatorGridLayouts[opr_count].addWidget(self.operatorTypeLabels[opr_count], 0, 0, 1, 1)
        self.operatorTypeLabels[opr_count].setText('Operator Type:')

        self.operatorPairStateLabels.append(QtWidgets.QLabel(self.scrollAreaWidgetContents_2))
        self.operatorPairStateLabels[opr_count].setObjectName(f'ooperatorPairStateLabels_{opr_count}')
        self.operatorGridLayouts[opr_count].addWidget(self.operatorPairStateLabels[opr_count], 1, 0, 1, 1)
        self.operatorPairStateLabels[opr_count].setText('Pair States:')

        self.operatorModelTypeLabels.append(QtWidgets.QLabel(self.scrollAreaWidgetContents_2))
        self.operatorModelTypeLabels[opr_count].setObjectName(f'ooperatorModelTypeLabels_{opr_count}')
        self.operatorGridLayouts[opr_count].addWidget(self.operatorModelTypeLabels[opr_count], 2, 0, 1, 1)
        self.operatorModelTypeLabels[opr_count].setText('Model Type:')

        self.operatorRadialLabelLabels.append(QtWidgets.QLabel(self.scrollAreaWidgetContents_2))
        self.operatorRadialLabelLabels[opr_count].setObjectName(f'ooperatorRadialLabelLabels_{opr_count}')
        self.operatorGridLayouts[opr_count].addWidget(self.operatorRadialLabelLabels[opr_count], 3, 0, 1, 1)
        self.operatorRadialLabelLabels[opr_count].setText('Label of radial parameters:')

        self.operatorRotCorrLabels.append(QtWidgets.QLabel(self.scrollAreaWidgetContents_2))
        self.operatorRotCorrLabels[opr_count].setObjectName(f'ooperatorRotCorrLabels{opr_count}')
        self.operatorGridLayouts[opr_count].addWidget(self.operatorRotCorrLabels[opr_count], 4, 0, 1, 1)
        self.operatorRotCorrLabels[opr_count].setText('Rot correction:')

        self.operatorTypeComboBox.append(QtWidgets.QComboBox(self.scrollAreaWidgetContents_2))
        self.operatorTypeComboBox[opr_count].setObjectName(f'ooperatorTypeComboBox_{opr_count}')
        self.operatorGridLayouts[opr_count].addWidget(self.operatorTypeComboBox[opr_count], 0, 1, 1, 2)
        self.operatorTypeComboBox[opr_count].addItems(self.opr_type_names)
        self.operatorTypeComboBox[opr_count].currentTextChanged.connect(
            lambda: self._on_operator_type_change(opr_label-1)
        )
        
        self.pairStates1SpinBox.append(QtWidgets.QSpinBox(self.scrollAreaWidgetContents_2))
        self.pairStates1SpinBox[opr_count].setObjectName(f'pairStates1SpinBox_{opr_count}')
        self.pairStates1SpinBox[opr_count].setMinimum(1)
        self.operatorGridLayouts[opr_count].addWidget(self.pairStates1SpinBox[opr_count], 1, 1, 1, 1)
        self.pairStates1SpinBox[opr_count].valueChanged.connect(
            lambda: self._on_operator_pair_state1_value_change(opr_label-1)
        )

        self.pairStates2SpinBox.append(QtWidgets.QSpinBox(self.scrollAreaWidgetContents_2))
        self.pairStates2SpinBox[opr_count].setObjectName(f'pairStates2SpinBox_{opr_count}')
        self.pairStates2SpinBox[opr_count].setMinimum(1)
        self.operatorGridLayouts[opr_count].addWidget(self.pairStates2SpinBox[opr_count], 1, 2, 1, 1)
        self.pairStates2SpinBox[opr_count].valueChanged.connect(
            lambda: self._on_operator_pair_state2_value_change(opr_label-1)
        )

        self.modelTypeComboBox.append(QtWidgets.QComboBox(self.scrollAreaWidgetContents_2))
        self.modelTypeComboBox[opr_count].setObjectName(f'modelTypeComboBox_{opr_count}')
        self.operatorGridLayouts[opr_count].addWidget(self.modelTypeComboBox[opr_count], 2, 1, 1, 2)
        self.modelTypeComboBox[opr_count].addItems(self.model_types)
        self.modelTypeComboBox[opr_count].currentTextChanged.connect(
            lambda: self._on_operator_model_change(opr_label-1)
        )

        self.labelRadParamsLineEdit.append(QtWidgets.QLineEdit(self.scrollAreaWidgetContents_2))
        self.labelRadParamsLineEdit[opr_count].setObjectName(f'labelRadParamsLineEdit_{opr_count}')
        self.operatorGridLayouts[opr_count].addWidget(self.labelRadParamsLineEdit[opr_count], 3, 1, 1, 2)
        self.labelRadParamsLineEdit[opr_count].editingFinished.connect(
            lambda: self._get_operator_label(opr_label-1)
        )

        self.labelRotCorrLineEdit.append(QtWidgets.QLineEdit(self.scrollAreaWidgetContents_2))
        self.labelRotCorrLineEdit[opr_count].setObjectName(f'labelRotCorrLineEdit_{opr_count}')
        self.operatorGridLayouts[opr_count].addWidget(self.labelRotCorrLineEdit[opr_count], 4, 1, 1, 2)
        self.labelRotCorrLineEdit[opr_count].setPlaceholderText('0.0')
        self.labelRotCorrLineEdit[opr_count].editingFinished.connect(
            lambda: self._get_operator_rot_corr(opr_label-1)
        )

        self.showRadParamsPushButton.append(QtWidgets.QPushButton(self.scrollAreaWidgetContents_2))
        self.showRadParamsPushButton[opr_count].setObjectName(f'showRadParamsPushButton_{opr_count}')
        self.showRadParamsPushButton[opr_count].setText('Show Radial Parameters')
        self.operatorGridLayouts[opr_count].addWidget(self.showRadParamsPushButton[opr_count], 5, 0, 1, 1)

        self.plotRadParamsPushButton.append(QtWidgets.QPushButton(self.scrollAreaWidgetContents_2))
        self.plotRadParamsPushButton[opr_count].setObjectName(f'plotRadParamsPushButton_{opr_count}')
        self.plotRadParamsPushButton[opr_count].setText('Plot Radial Function')
        self.operatorGridLayouts[opr_count].addWidget(self.plotRadParamsPushButton[opr_count], 5, 1, 1, 2)

        self.gridLayout_8.addLayout(self.operatorGridLayouts[opr_count], opr_count, 0, 1, 1)
        self.operatorTypeLabels[opr_count].show()

def main():

    app = QtWidgets.QApplication(sys.argv)
    app.setStyle('Fusion')

    window = Diatom()

    screen_resolution = app.desktop().screenGeometry()
    width, height = screen_resolution.width(), screen_resolution.height()

    if height < 920:
        window.setMinimumWidth(width*0.6)
        window.setMinimumHeight(height*0.6)
        window.resize(width*0.8, height*0.8)
    else:
        pass

    window.show()

    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
