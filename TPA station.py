import sys
from time import sleep
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
from PyQt5.QtWidgets import QApplication, QWidget, QListWidgetItem, QFileDialog, QMessageBox
from PyQt5 import uic
from PyQt5.QtGui import QPixmap
from Ressources_scripts.Data_processing_tpa_calculation import Read_Data_File as rdf
from Ressources_scripts.Data_processing import get_all_samples, save, get_samples
from pathlib import Path

class MplCanvas(FigureCanvas):
    def __init__(self, parent=None):
        fig = Figure()
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)
        fig.tight_layout()

class MyApp(QWidget):

    def __init__(self):
        super().__init__()
        uic.loadUi('TPA station.ui', self)
        self.setWindowTitle('TPA station')
        self.processed_data = {}
        self.label_Logo_NMP.setPixmap(QPixmap("Logo NMP.png"))

        """
        This app contains:
        10 buttons:
            *btn_get_files_2: ask the user to give data folder path
            *btn_set_sample: do nothing
            *btn_set_sample_wavelength: select sample & wavelength then show data
            *btn_simulation: remove checked data and recalculate curves fitting and 2PA, DATA ARE NOT SAVE
            *btn_remove_data: remove checked data and recalculate curves fitting and 2PA, CHANGE THE DATA SOURCE FILES, CANNOT BE UNDONE
            *btn_show_sigma_2: display 2PA spectrum if necessary of the current dye
            *btn_figure_2p_vs_1p: save the emission spectra of 2PA and 1PA as a .png file far current dye & wavelength         
            *btn_refresh: refresh sample & wavelength
            *btn_change_sample_data: allow user to change concentration, phi, solvent and emission spectrum
            *btn_3d_plot: plot (wavelength, fluorescence, power) plot of selected sample
        34 labels:
            *10 associated to fluorescence: label_f_1 to label_f_8 and label_f_18 & label_f_17 (for 9th & 10th position respectively)
            *10 associated to laser_power: label_P_1 to label_P_10
            *10 associated to fluorescence/(laser_power)²: label_Quad_1 to label_Quad_10
            *4 associated to variance and standard deviation
        4 plotwidgets:
            *S2Spectrum: show two-photon absorption (2PA) spectrum of the current dye
            *FPSpectrum: show fluorescence = f(laser power) with curve fitting and equation of the current dye & wavelength
            *FluoIntensitySpectrum: show all fluorescence spectra of the current dye & wavelength at different laser power
            *FPSpectrum_2: show log(fluorescence) = f(log(laser_power²)) with curve fitting and slope value of the current dye & wavelength
        20 qcheckboxes:
            *10 associated to simulation: chk_simulation_1 to chk_simulation_10
            *10 associated to remove: chk_remove_1 to chk_remove_10
        2 qlists:
            *ls_list_of_all_samples: show all measured samples and allow user to select one of them
            *ls_done_wavelengths: show all measured wavelength for the selected sample and allow user to select one of them
        """

        """Variables"""
        self.user_info_measurement = []
        self.order_of_measure = []
        self.sigma2_dict = {}
        self.wavelength_list = []
        self.sigma2_list = []
        self.sigma2phi_list = []
        self.experimental_data = {}
        self.counter_results = {}
        self.root_dict = {}
        self.done_wavelengths = []

        """Button"""
        self.btn_get_files_2.clicked.connect(self.read_data_from_folder)
        self.btn_remove_data.clicked.connect(self.remove_data)
        self.btn_refresh.clicked.connect(self.refresh_data)
        self.btn_set_sample_wavelength.clicked.connect(self.plot_s2_from_signal)
        self.btn_show_sigma_2.clicked.connect(self.plot_s2_from_signal)
        self.btn_simulation.clicked.connect(self.simulation)
        self.btn_figure_2p_vs_1p.clicked.connect(self.generate_figure_2p_vs_1p)
        self.btn_change_sample_data.clicked.connect(self.change_sample_data)
        self.btn_3d_plot.clicked.connect(self.plot_3d)

    def read_data_for_3d(self, sample):
        x=[]
        y=[]
        z=[]
        print("Read_data_for_3d")

        for wavelength in self.list_of_wavelength:
            print("|"*1000)
            print(wavelength)
            try :
                fluo, power = save()._get_quick_f_and_pp(sample=sample, wavelength=wavelength, root= self.folder_root, number_of_measure=self.power_pitch)
                for n in range(len(fluo)):
                    y.append(power[n])
                    z.append(fluo[n])
                    x.append(int(wavelength))
                print(len(x), len(y), len(z))
                print(x, y, z)
            except:
                print("Error in wavelength reading")

        x = np.array(x)
        y=np.array(y)
        z=np.array(z)
        return x, y, z

    def plot_3d(self):
        sample= self.select_sample()
        print(sample)
        x, y, z = self.read_data_for_3d(sample=sample)
        print("3d")
        fig, ax = plt.subplots()
        #ax = plt.axes(projection='3d')
        """ax.set_xlabel("Wavelength / nm")
        ax.set_zlabel("F / cps")
        ax.set_ylabel("P² / mA")"""
        y_2d, z_2d=np.meshgrid(y, z)
        #ax.scatter(x, y, z, c=z, cmap='BrBG', linewidth=1)

        #print(x_2d.shape, y_2d.shape, z_2d.shape)
        #ax.plot_wireframe(x, y, z_2d, rstride=10, cstride=10)

        #ax.plot_trisurf(x, y, z, cmap='viridis')
        #ax.plot_surface(x, y, z_2d, cmap='jet')
        ax.pcolormesh(x, y, z_2d)
        #fig.savefig(DATA / f"{name}.png")
        #plt.draw()
        """for angle in range(-150, 0):
            ax.view_init(elev=10, azim=angle, vertical_axis="z")
            plt.draw()
            plt.pause(.01)"""
        plt.show()


    def change_sample_data(self):
        sample= self.select_sample()
        dlg = QMessageBox(self)
        dlg.setWindowTitle(f"{sample}")
        dlg.setText("What do you want to change?")
        buttonoptionE = dlg.addButton("Cancel", QMessageBox.ButtonRole.AcceptRole)
        buttonoptionD = dlg.addButton("Emission spectrum", QMessageBox.ButtonRole.AcceptRole)
        buttonoptionC = dlg.addButton("Solvent", QMessageBox.ButtonRole.AcceptRole)
        buttonoptionB = dlg.addButton("Phi", QMessageBox.ButtonRole.AcceptRole)
        buttonoptionA = dlg.addButton("Concentration", QMessageBox.ButtonRole.AcceptRole)
        original_data= save().read_samples_informations(self.folder_root)
        dlg.exec()

        new_data=0#temporaire à enlever

        if dlg.clickedButton() == buttonoptionA:
            """
            Ici faire :
            1-poper une boîte de dialog avec une entrée user => ici vérifier que l'user rentre un float
            2-Sauver dans une variable
            3-updater le dico original_data avec la new_data=> rôle de update_samples_informations
            """

            QMessageBox.information(self, "Information", f"Concentration of {sample} changed")
            save().update_samples_informations(original_data=original_data, new_data=new_data, type_of_data="Concentration")
            print("option A")
        elif dlg.clickedButton() == buttonoptionB:
            QMessageBox.information(self, "Information", f"Phi of {sample} changed")
            save().update_samples_informations(original_data=original_data, new_data=new_data,
                                               type_of_data="Phi")
            print("option B")
        elif dlg.clickedButton() == buttonoptionC:
            QMessageBox.information(self, "Information", f"Solvent of {sample} changed")
            save().update_samples_informations(original_data=original_data, new_data=new_data,
                                               type_of_data="Solvent")
            print("option C")
        elif dlg.clickedButton() == buttonoptionD:
            QMessageBox.information(self, "Information", f"Emission spectrum of {sample} changed")
            save().update_samples_informations(original_data=original_data, new_data=new_data,
                                               type_of_data="Emission spectrum")
            print("option D")

        self.show_data_of_sample(simulation=False)

    def read_data_from_folder(self):
        temp = QFileDialog.getExistingDirectory(caption='Please select the folder Datas')
        self.folder_root = Path(temp)
        self.sample_root = self.folder_root / "Metadata" / "samples.json"
        self.sample_data_root = self.folder_root / "Metadata" /"sample_data.json"
        print(self.folder_root)
        print(self.sample_data_root)
        self.refresh_data()
        self.get_files()

    def generate_figure_2p_vs_1p(self):
        wavelength= self.select_wavelength()
        sample= self.select_sample()
        n, data = save().extract_data_from_experimental(sample=sample, wavelength=wavelength, root=self.folder_root)
        for i in range(len(n)):
            x= data[f"{wavelength}"][f"{i}"]["wavelength"]
            y= data[f"{wavelength}"][f"{i}"]["intensity - dark"]
            intensity_max=max(y)
            print(intensity_max)
            temp_y=[]
            temp_y2=[]
            for z in y:
                z /= intensity_max
                "Normalization step"
                temp_y.append(z)
            x2, y2= save().read_column_from_excel_files(sample=sample, root=self.folder_root)
            intensity_max_2=max(y2)
            for h in y2:
                h /= intensity_max_2
                temp_y2.append(h)
                "Normalization step"
            y=np.array(temp_y)
            y2=np.array(temp_y2)
            print(max(y))
            print(max(y2))
            save().figure_emission_2p_vs_1p(x1=x, y1=y, x2=x2, y2=y2, sample=sample, wavelength=f"{wavelength}_nm", root=self.folder_root, power=i, wavenumber=False)
            temp_y.clear()
            temp_y2.clear()
            for z in x:
                z = 1/(z*1e-7)
                temp_y.append(z)
            for h in x2:
                h = 1/(h*1e-7)
                temp_y2.append(h)
            x= np.array(temp_y)
            x2=np.array(temp_y2)
            save().figure_emission_2p_vs_1p(x1=x, y1=y, x2=x2, y2=y2, sample=sample, wavelength=f"{wavelength}_cm-1", root=self.folder_root, power=i, wavenumber= True)


    def get_files(self):
        self.sample_info = save().read_samples_informations(root=self.folder_root)
        print("def GET_FILES")
        print("self.sample_info: ", self.sample_info)
        samples= get_samples(root= self.sample_root)
        print("sample: ", samples)
        self.list_of_sample_data ={}
        self.list_of_wavelength=[]
        self.list_of_sample_data[samples[0]] = save().read_process(samples[0], simulation=False, root=self.folder_root)
        for wave in self.list_of_sample_data[samples[0]]:
            self.list_of_wavelength.append(int(wave))
            print("self.list_of_wavelength: ", self.list_of_wavelength, type(self.list_of_wavelength))
        self.refresh_wavelength(self.list_of_wavelength)
        self.power_pitch = self.sample_info["i"]
        print("i= ", self.power_pitch)
        self.fluorescence_ref_dye = self.sample_info["TPEF"]["name"]
        print("~"*1000)
        print(type(self.fluorescence_ref_dye))

        for sample in samples:
            print(sample)
            self.processed_data.update({f"{sample}":save().read_process(sample= sample, simulation=False, root=self.folder_root)})
            save().copy_processed_data_for_simulation(sample, root=self.folder_root)

    def simulation(self):
        list_of_measure_to_remove = []
        for selected_data in self.ls_measure_data.selectedItems():
            a = self.ls_measure_data.row(selected_data)


            list_of_measure_to_remove.append(a)
        for selected_sample in self.ls_list_of_all_samples.selectedItems():
            sample = selected_sample.data(0x0100)


        for selected_wavelength in self.ls_done_wavelengths.selectedItems():
            wavelength = selected_wavelength.data(0x0100)

        save().copy_data_for_simulation(sample=sample, wavelength= wavelength, root=self.folder_root)
        save().copy_processed_data_for_simulation(sample=sample, root=self.folder_root)
        for i in list_of_measure_to_remove:
            print(list_of_measure_to_remove)
            print("simulate without : ", i, sample, wavelength)
            save().remove_data(sample=sample, wavelength=wavelength, measure=i, simulation=True, root=self.folder_root)

        self.tpef_calculation_for_single_measure(sample=sample, wavelength=wavelength, simulation=True)
        self.plot_s2_from_signal(simulation=True)

    def graph_data(self, x, y,x2, y2, view, value, order_process, x_log, y_log, x_log2, y_log2):
        x = np.array(x)
        y = np.array(y)
        if view == "fluo":
            self.FluoIntensitySpectrum.plot(x, y, show= True)
        elif view == "S2F":
            self.S2Spectrum.clear()
            self.S2Spectrum.plot(x, y, show = True)
        elif view == "quad":
            self.FPSpectrum.clear()
            self.FPSpectrum.addLegend()
            self.FPSpectrum.setLabel('left', text= "F / cps")
            self.FPSpectrum.setLabel('bottom', text="P²")
            self.FPSpectrum.plot(x, y, pen= "k", symbol= "o", symbolPen="g",symbolBrush=0.5, show= True)
            self.FPSpectrum.plot(x2, y2, show=True)
            self.FPSpectrum.setLabel('top', text=f"y= {value}")
            self.FPSpectrum_2.clear()
            self.FPSpectrum_2.addLegend()
            self.FPSpectrum_2.setLabel('left', text= "log(F) / cps")
            self.FPSpectrum_2.setLabel('bottom', text="log(P²)")
            self.FPSpectrum_2.plot(x_log, y_log, pen= "k", symbol= "o", symbolPen="g",symbolBrush=0.5, show= True)
            self.FPSpectrum_2.plot(x_log2, y_log2, show=True)
            self.FPSpectrum_2.setLabel('top', text=f"Slope = {order_process} (ref: 1.0)")

    def refresh_data(self):
        self.ls_list_of_all_samples.clear()
        samples = get_samples(self.sample_root)
        for sample in samples:
            lw_item = QListWidgetItem(sample.sample_name)
            lw_item.setData(0x0100, sample)  # 0x0100 => QtCore.Qt.UserRole
            self.ls_list_of_all_samples.addItem(lw_item)
            print(sample)

    def refresh_wavelength(self, list_of_wavelength):
        self.ls_done_wavelengths.clear()
        for wavelength in list_of_wavelength:
            wavelength= int(wavelength)
            lw_item = QListWidgetItem(f"{wavelength}")
            lw_item.setData(0x0100, wavelength)  # 0x0100 => QtCore.Qt.UserRole
            self.ls_done_wavelengths.addItem(lw_item)
            print(wavelength)

    def select_wavelength(self):
        for selected_wavelength in self.ls_done_wavelengths.selectedItems():
            wavelength = selected_wavelength.data(0x0100)
            print("The selected wavelength is: ", wavelength)
        return wavelength

    def select_sample(self):
        for selected_sample in self.ls_list_of_all_samples.selectedItems():
            sample = selected_sample.data(0x0100)
            print("The selected sample is : ", sample)
            return sample

    def plot_s2_from_signal(self, simulation):
        self.S2Spectrum.clear()
        self.show_data_of_sample(simulation=simulation)
        sample= self.select_sample()
        data= save().read_process(sample= sample, simulation=simulation, root=self.folder_root)
        sigma_2= []
        tpef_wavelength= []
        for wavelength in data:
            sigma_2.append(int(data[f"{wavelength}"]["S2"]))
            print(type(sigma_2))
            tpef_wavelength.append(int(wavelength))
            print(type(wavelength))
        x = np.array(tpef_wavelength)
        y= np.array(sigma_2)
        print("Wavelength:", x)
        print("S2:", y)
        self.S2Spectrum.setLabel('left', text="Sigma 2 / GM")
        self.S2Spectrum.setLabel('bottom', text="Wavelength / nm")
        self.S2Spectrum.plot(x, y, pen="g", symbol="o", symbolPen="m", symbolBrush= 0.5, name= f"{sample}", show=True)
        try:
            save().figure(x=x, y=y, name=f"TPEF spectra of {sample}", sample=sample, root=self.folder_root)
        except ImportError:
            print(f"Import Error for the figure TPEF spectra of {sample}")

    def save_quad_from_signal(self, square_power, fluo, fit, sample, wavelength):
        save().figure_quadraticity(x1=square_power, y1=fluo, x2=square_power, y2=fit,
                                   sample=sample, wavelength=wavelength)

    def save_fluo_vs_power(self, i, temporary_wavelength_list, temporary_intensity_list, sample, wavelength):
        save().figure_intensity_vs_power(number_of_measure=i, x=temporary_wavelength_list,
                                         y=temporary_intensity_list, sample=sample, wavelength=wavelength)

    def clear_plot_from_signal(self):
        print("#"*100)
        print("Clear graphs")
        print("#" * 100)
        self.S2Spectrum.clear()
        self.FluoIntensitySpectrum.clear()
        self.FPSpectrum.clear()
        self.S2Spectrum.addLegend()
        self.FluoIntensitySpectrum.addLegend()
        self.FPSpectrum.addLegend()

    def plot_from_signal(self, x, y, square_power_dark_corr, full_corrected_area, i, power):

        dictionary_color = {"0": "w", "1": "b", "2": "c", "3": "g", "4": "y", "5": "m", "6": "r", "7": "w",
                            "8": "b", "9": "c", "10":"g"}

        try:
            power= int(power)
        except ValueError:

            print("error in plot")

        x = np.array(x)
        y = np.array(y)
        self.FluoIntensitySpectrum.addLegend()
        self.FluoIntensitySpectrum.setLabel('left', text="Fluorescence intensity / cps")
        self.FluoIntensitySpectrum.setLabel('bottom', text="Wavelength / nm")
        if i > 10:

            data=str(i)
            data_temp=data[1:]
            i=int(data_temp)

        self.FluoIntensitySpectrum.plot(x, y, name= f"{power} \n", pen= dictionary_color[f"{i}"], show=True)

    def save_plot_from_signal(self, x, y, wavelength, sample):
        save().figure(x=x, y=y, name=f"Dark at {wavelength}", sample=sample)

    def update_dictionary_processed_data(self, sample, wavelength, slope, coeff, log_square_fluo, log_square_power, fit,
                                         full_corrected_area, sigma2_phi, sigma2, list_of_fl, list_of_pp):
        self.processed_data.update({
            sample:
                {
                    f"{wavelength}":
                        {
                            "quadraticity":
                                {
                                    'slope': slope,
                                    'coeff': coeff,
                                    'log(F)': log_square_fluo,
                                    'log(P²)': log_square_power,
                                    'fit': fit,
                                    'F': list_of_fl,
                                    'PP': list_of_pp
                                },
                            'full_corrected_area': full_corrected_area,
                            "S2F": sigma2_phi,
                            "S2": sigma2
                        }
                }
        })

    def tpef_calculation_for_single_measure(self, sample, wavelength, simulation):
        list_of_fully_corrected_area= []
        list_of_square_dark_corr_power = []
        data = save()._get_data_for_new_mesure(sample= sample, wavelength= wavelength, simulation=simulation, root=self.folder_root)
        for key_wavelength in data:
            print("Starting TPEF calulation.")
            print(type(data))
            for key_measure in data[f"{key_wavelength}"]:
                if key_measure != "dark":
                    print(key_measure)
                    print(type(key_measure))
                    print(data[f"{key_wavelength}"].keys())
                    print(data[f"{key_wavelength}"][f"{key_measure}"].keys())
                    print(data[f"{key_wavelength}"][f"{key_measure}"]["raw power - dark"])
                    print(data[f"{key_wavelength}"][f"{key_measure}"]["fully corrected area"])
                    power = data[f"{key_wavelength}"][f"{key_measure}"]["raw power - dark"]
                    print("power = ", power)
                    square_power = power*power
                    print('Square power = ', square_power)
                    list_of_square_dark_corr_power.append(square_power)
                    list_of_fully_corrected_area.append(data[f"{key_wavelength}"][f"{key_measure}"]["fully corrected area"])

        area, correlation_coeff, slope, coeff, log_square_fluo, log_square_power, fit = rdf()._Quad_Log(list_of_fully_corrected_area,
                                                                               list_of_square_dark_corr_power)

        print("+" * 100)
        print(self.fluorescence_ref_dye, type(self.fluorescence_ref_dye))
        print(sample, type(sample))
        print("+"*100)
        sigma2_phi, sigma2 = rdf()._Calcul_TPA(sample_info=self.sample_info,
                                               processed_data=self.processed_data,
                                               sample_name=str(sample),
                                               reference_name=self.fluorescence_ref_dye,
                                               wavelength=int(wavelength),
                                               simulation=simulation,
                                               root=self.folder_root)
        print("Update dictionnary")
        self.update_dictionary_processed_data(sample=sample, wavelength=wavelength, slope=slope,
                                              coeff=coeff, log_square_fluo=log_square_fluo,
                                              log_square_power=log_square_power, fit=fit,
                                              full_corrected_area=area,
                                              sigma2_phi=sigma2_phi, sigma2=sigma2, list_of_fl=list_of_fully_corrected_area, list_of_pp=list_of_square_dark_corr_power)
        print("self.processed_data: ", self.processed_data)
        sleep(0.1)

        file_process = save().save_process(dictionary=self.processed_data[sample], sample=sample, state=False, simulation=simulation, root=self.folder_root)
        self.root_dict.update({sample: file_process})
        if simulation == False:
            print("simulation is false! ")
            file_exp = save()._get_file_for_new_mesure(sample= sample, wavelength= wavelength, root= self.folder_root)
            #save().figure_quadraticity(x1=log_square_power, y1=log_square_fluo, x2=log_square_power, y2=fit, sample=sample, wavelength=wavelength, root=self.folder_root)
            save().save_experimental_to_excel(file=file_exp, sample=sample, root= self.folder_root)
            save().save_process_data_to_excel(file=file_process, sample=sample, root= self.folder_root)
            self.root_dict.update({sample: file_process})

            file = self.root_dict[sample]
            data_for_sigma_2 = save().read_json_data(file)
            print(20 * "#")
            print(sample)
            print(data_for_sigma_2.keys())
            print(data_for_sigma_2[f"{wavelength}"].keys())
            print("S2F = ", data_for_sigma_2[f"{wavelength}"]["S2F"])
            x = []
            y = []
            for wavelength in self.list_of_wavelength: # self.wavelength_list:
                print(data_for_sigma_2.keys())
                print(data_for_sigma_2.values)
                print(data_for_sigma_2[f"{wavelength}"])
                x.append(wavelength)
                y.append(data_for_sigma_2[f"{wavelength}"]["S2F"])
            try :
                save().figure(x=x, y=y, name=f"TPEF spectra of {sample}", sample=sample, root= self.folder_root)
            except ImportError:
                print(f"Import Error for the figure TPEF spectra of {sample}")

    def remove_data(self):
        list_of_measure_to_remove = []
        for selected_data in self.ls_measure_data.selectedItems():
            a= self.ls_measure_data.row(selected_data)

            list_of_measure_to_remove.append(int(a))

        for selected_sample in self.ls_list_of_all_samples.selectedItems():
            sample = selected_sample.data(0x0100)
        for selected_wavelength in self.ls_done_wavelengths.selectedItems():
            wavelength = selected_wavelength.data(0x0100)
        for i in list_of_measure_to_remove:
            save().remove_data(sample=sample, wavelength=wavelength, measure=i, simulation=False,root=self.folder_root)
        self.tpef_calculation_for_single_measure(sample=sample, wavelength=wavelength, simulation= False)
        self.plot_s2_from_signal(simulation=True)
        print("** Remove data point successful ! **")

    def show_sample_name(self):
        samples = get_all_samples()
        for sample in samples:
            lw_item = QListWidgetItem(sample.sample_name)
            lw_item.setData(0x0100, sample)  # 0x0100 => QtCore.Qt.UserRole
            self.ls_list_of_samples.addItem(lw_item)


    def show_data_of_sample(self, simulation):
        self.ls_measure_data.clear()
        number_of_measure = self.power_pitch
        for selected_sample in self.ls_list_of_all_samples.selectedItems():
            sample = selected_sample.data(0x0100)
        for selected_wavelength in self.ls_done_wavelengths.selectedItems():
            wavelength = selected_wavelength.data(0x0100)
            if simulation == False:
                file_root = self.folder_root / f"{sample}" / f"experimental_data {sample}_{wavelength}.json"
                list_of_fl, list_of_pp = save()._get_quick_f_and_pp(sample=sample, wavelength=wavelength,
                                                                    root=self.folder_root, number_of_measure= self.power_pitch)
                data = save().read_process(sample=sample, simulation=simulation, root=self.folder_root)

                slope= data[f"{wavelength}"]["quadraticity"]['slope']
                coeff=data[f"{wavelength}"]["quadraticity"]['coeff']
                log_square_fluo=data[f"{wavelength}"]["quadraticity"]["log(F)"]
                log_square_power=data[f"{wavelength}"]["quadraticity"]['log(P²)']
                fit=data[f"{wavelength}"]["quadraticity"]['fit']
                full_corrected_area=data[f"{wavelength}"]['full_corrected_area']
                sigma2_phi=data[f"{wavelength}"]['S2F']
                sigma2=data[f"{wavelength}"]['S2']


                data.update(
                        {
                            f"{wavelength}":
                                {
                                    "quadraticity":
                                        {
                                            'slope': slope,
                                            'coeff': coeff,
                                            'log(F)': log_square_fluo,
                                            'log(P²)': log_square_power,
                                            'fit': fit,
                                            'F': list_of_fl,
                                            'PP': list_of_pp
                                        },
                                    'full_corrected_area': full_corrected_area,
                                    "S2F": sigma2_phi,
                                    "S2": sigma2
                                }
                        }
                )

                save().save_process(dictionary=data, sample=sample, state=False,
                                               simulation=simulation, root=self.folder_root)
            elif simulation == True:
                file_root = self.folder_root / f"{sample}" / f"experimental_data {sample}_{wavelength}_temp.json"
            data = save().read_json_data(file_root)
            print(file_root)

            temp_quadraticity_data=[]

            for i in range(number_of_measure):
                try:
                    print(data.keys())
                    print(data[f"{wavelength}"].keys())
                    fluorescence_intensity = data[f"{wavelength}"][f"{i}"]['fully corrected area']
                    power = data[f"{wavelength}"][f"{i}"]['raw power - dark']
                    square_power = power*power
                    quadraticity = (float(fluorescence_intensity)/float(square_power))
                    square_power = round(square_power, 5)

                    temp_quadraticity_data.append(quadraticity)
                    self.ls_measure_data.addItem(f"{i} \t {fluorescence_intensity:.2E} \t {square_power} \t {quadraticity:.2E}")
                except KeyError:
                    print("key error _show data sample !")

        temp_quadraticity_data_array=np.array(temp_quadraticity_data)

        mean_data=np.average(a=temp_quadraticity_data)
        print("Mean = ", mean_data)
        self.Qlabel_mean.setText(f"{mean_data:.1E}")
        self.QLabel_mean_text.setText("+/-")

        variance_data = np.var(a=temp_quadraticity_data_array)
        print("Variance = ", variance_data)
        self.Qlabel_variance.setText(f"{variance_data:.1E}")

        standard_deviation_data=np.std(a=temp_quadraticity_data_array)
        print("Standard deviation = ", standard_deviation_data)
        self.Qlabel_standard_deviation.setText(f"{standard_deviation_data:.1E}")

        s2f, x_wavelength, log_f, log_pp, list_fl, list_pp = save().extract_process_data(sample=sample,
                                                                       wavelength=wavelength, simulation=simulation, root=self.folder_root)

        x2, y2, value, order_process, x_log, y_log, x_log2, y_log2= save().fit_curve(x= np.array(list_pp), y=np.array(list_fl))

        self.graph_data(x=list_pp, y=list_fl, x2=x2, y2= y2, value= value, order_process= order_process, view="quad", x_log=x_log, y_log=y_log, x_log2=x_log2, y_log2=y_log2)

        self.FluoIntensitySpectrum.clear()
        i=0
        data_bis= save().extract_data_from_experimental_bis(sample= sample, wavelength= wavelength, root= self.folder_root)
        for key in data_bis:
                for power in data_bis[f"{key}"]:
                    if power != "dark":
                        i += 1
                        print("*"*100)
                        print(power)
                        x, y = save()._get_intensity_vs_wavelength_data(data_bis[f"{key}"][f"{power}"])
                        print(len(x))
                        print(len(y))
                        print("*" * 100)
                        self.plot_from_signal(x= x, y=y, square_power_dark_corr="", full_corrected_area="", i= i, power= i)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyleSheet('''
        QWidget {
        font-size: 30 px;
        }
    ''')

    myApp = MyApp()
    myApp.show()
    sys.exit(app.exec())

    try:
        sys.exit(app.exec())
    except SystemExit:
        print('Closing Windows')