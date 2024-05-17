import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QTreeView, QPlainTextEdit, QCheckBox
from PyQt5.QtWidgets import QPushButton, QVBoxLayout, QWidget, QFileSystemModel, QLabel, QHBoxLayout
from PyQt5.QtWidgets import QProgressBar, QHeaderView, QLineEdit, QGroupBox
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import QDir, Qt
from PyQt5 import QtCore
from Log import Log

from Fouille import Fouille
from Tree import Root
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

from PyQt5.QtCore import QThread, pyqtSignal
import os

class ProcessRecordsThread(QThread):
    progressSignal = pyqtSignal(int, int, int, int)

    def __init__(self, fouille, arbre, selectedRegions):
        super().__init__()
        self.fouille = fouille
        self.arbre = arbre
        self.selectedRegions = selectedRegions
        self.percentage = 0
        self.nc_total = 0
        self.processed_ncs = 0
        self.total_organism = 0
        self.processed_organism = 0

    def run(self):
        self.fouille = Fouille(self.selectedRegions)
        cache_name = "logs/aborted-parsing.txt"
        abort_log = Log(cache_name)
        self.fouille.peuplement(self.arbre, abort_log, [ l[:-1] for l in abort_log.gen_list() ],  self.updateProgress)
        self.progressSignal.emit(
            self.fouille.get_total_organisms(),
            self.fouille.get_processed_organisms(),
            self.fouille.get_total_ncs(),
            self.fouille.get_processed_ncs()
        )


    def updateProgress(self, total_organism, processed_organism,  nc_total, processed_ncs):
        self.progressSignal.emit(total_organism, processed_organism, nc_total, processed_ncs)
        self.total_organism = total_organism
        self.processed_organism = processed_organism
        self.nc_total = nc_total
        self.processed_ncs = processed_ncs

    
    def getTotalOrganism(self):
        return self.total_organism
    
    def getProcessedOrganism(self):
        return self.processed_organism

    def getNcTotal(self):
        return self.nc_total
    
    def getProcessedNcs(self):
        return self.processed_ncs

class CustomFileSystemModel(QFileSystemModel):
    def __init__(self):
        super().__init__()

    def data(self, index, role):
        if role == QtCore.Qt.DecorationRole and index.isValid():
            fileInfo = self.fileInfo(index)
            if fileInfo.isDir():
                if index.column() == 0:
                    return QIcon('./img/pink-folder.png')  
            else:
                if index.column() == 0:
                    return QIcon('./img/pink-file.png') 
        return super().data(index, role)

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Record Processing")
        self.setGeometry(100, 100, 1600, 1000)
        self.setStyleSheet("background-color: white;") 

        self.arbre = Root()
        self.fouille = Fouille([])
        self.percentage = 0
        self.processRecordsThread = None
        self.organismPercentage = 0
        self.total_ncs = 0 # for current subtree

        self.setupUi()

    def setupUi(self):
        self.treeView = QTreeView()
        self.treeView.setRootIsDecorated(True)
        self.treeView.setAlternatingRowColors(True)
        self.model = QFileSystemModel()
        self.model.setRootPath(QDir.rootPath())
        self.treeView.setModel(self.model)
        self.treeView.header().hideSection(1)
        self.treeView.header().hideSection(2)
        self.treeView.header().hideSection(3)
        self.treeView.header().setSectionResizeMode(0, QHeaderView.Stretch)
        self.treeView.setStyleSheet(            
            "QTreeView {"
            "border: 2px solid pink;"
            "border-radius: 5px;"
            "background-color: white;"
            "}")

        self.textEdit = QPlainTextEdit()
        self.textEdit.setReadOnly(True)
        self.processButton = QPushButton("Process Records")
        self.logsTextEdit = QPlainTextEdit()
        self.logsTextEdit.setReadOnly(True)
        self.organismProgressLabel = QLabel(f"Ogranisms processed: {self.organismPercentage:.3f}%")
        self.progressLabel = QLabel(f"NCs processed: {self.percentage:.3f}%")
        self.progressBar = QProgressBar(self)
        self.progressBar.setMaximum(10000)
        self.progressBar.setFormat("{:.3f}%".format(self.percentage))
        self.regionInput = QLineEdit()
        groupBox = QGroupBox()
        self.setStyleSheet("QLabel { background-color: pink; border-radius: 5px; padding: 5px; }")
        groupBox.setStyleSheet(
            "QGroupBox {"
            "border: 2px solid pink;"
            "border-radius: 5px;"
            "background-color: white;"
            "}"
        )
        leftLayout = QVBoxLayout()
        leftLayout.addWidget(self.treeView)

        headerLayout = QVBoxLayout()
        headerLayout.addWidget(QLabel("<b>Select Regions:</b>"))
        self.regionCheckBoxes = []
        regions = ['CDS', 'intron', 'telomere', 'centromere', 'mobile_element', 'ncRNA', 'rRNA', 'tRNA', '3\'UTR', '5\'UTR']
        for region in regions:
            checkBox = QCheckBox(region)
            headerLayout.addWidget(checkBox)
            self.regionCheckBoxes.append(checkBox)
        headerLayout.addWidget(QLabel("<b>Type a region:</b>"))
        headerLayout.addWidget(self.regionInput)
        groupBox.setLayout(headerLayout)
        rightLayout = QVBoxLayout()
        rightLayout.addWidget(groupBox) 

        groupBoxLogs = QGroupBox()
        groupBoxLogs.setStyleSheet(
            "QGroupBox {"
            "border: 2px solid pink;"
            "border-radius: 5px;"
            "background-color: white;"
            "}"
        )
        LogsLayout = QVBoxLayout()
        LogsLayout.addWidget(QLabel("<b>Logs:</b>"))
        LogsLayout.addWidget(self.logsTextEdit)
        groupBoxLogs.setLayout(LogsLayout)
        rightLayout.addWidget(groupBoxLogs)

        progressLayout = QVBoxLayout()
        progressLayout.addWidget(self.organismProgressLabel)
        progressLayout.addWidget(self.progressLabel)
        progressLayout.addWidget(self.progressBar)
        progressLayout.addWidget(self.processButton)
        groupBoxProgress = QGroupBox()
        groupBoxProgress.setLayout(progressLayout)
        groupBoxProgress.setStyleSheet(
            "QGroupBox {"
            "border: 2px solid pink;"
            "border-radius: 5px;"
            "background-color: white;"
            "}"
        )
        rightLayout.addWidget(groupBoxProgress)

        mainLayout = QHBoxLayout()
        mainLayout.addLayout(leftLayout)
        mainLayout.addLayout(rightLayout)
        mainLayout.addWidget(self.textEdit)
        mainLayout.setStretch(0, 2)
        mainLayout.setStretch(1, 1) 
        mainLayout.setStretch(2, 1)
        self.textEdit.setStyleSheet(            
            "QPlainTextEdit {"
            "border: 2px solid pink;"
            "border-radius: 5px;"
            "background-color: white;"
            "}")

        widget = QWidget()
        widget.setLayout(mainLayout)
        self.setCentralWidget(widget)

        # Connect signals
        self.treeView.clicked.connect(self.displayFileContent)
        self.processButton.clicked.connect(self.processRecords)

        self.populateTreeView()

    def populateTreeView(self):
        model = CustomFileSystemModel()
        root_path = os.path.join(QDir.currentPath(), "Results") # adapted for operating systems
        model.setRootPath(root_path)
        self.treeView.setModel(model)
        self.treeView.setRootIndex(model.index(root_path))

    def displayFileContent(self, index):
        model = self.treeView.model()
        filePath = model.filePath(index)
        if model.isDir(index):
            self.textEdit.clear()
            relative_path = os.path.relpath(filePath, os.path.join(QDir.currentPath(), "Results"))
            relative_path = "Results/"+relative_path
            print("Relative path:", relative_path)
            subtree = self.arbre.getSubtreeFromPath(relative_path) # Retrieve the subtree corresponding to the clicked item
            self.arbre.print_tree(subtree)
            return subtree
        else:
            with open(filePath, 'r') as file:
                content = file.read()
                self.textEdit.setPlainText(content)
            return None

    def processRecords(self):
        try:
            selectedRegions = []
            for checkBox in self.regionCheckBoxes:
                if checkBox.isChecked():
                    selectedRegions.append(checkBox.text())
            userRegion = self.regionInput.text()
            if userRegion:
                selectedRegions.append(userRegion)

            print("Selected Regions:", selectedRegions)
            # Call the function to process records with selected regions
            self.processRecordsThread = ProcessRecordsThread(self.fouille, self.arbre, selectedRegions)
            self.processRecordsThread.progressSignal.connect(self.updateProgress)
            self.processRecordsThread.start()
        except Exception as e:
            print("Error processing records:", e)

    def updateProgress(self, total_organism, processed_organism, nc_total, processed_ncs):
        total_organism = self.processRecordsThread.getTotalOrganism()
        processed_organism = self.processRecordsThread.getProcessedOrganism()
        nc_total = self.processRecordsThread.getNcTotal()
        processed_ncs = self.processRecordsThread.getProcessedNcs()
        
        self.organismPercentage = processed_organism/total_organism*100
        self.percentage = processed_ncs/nc_total*100
        self.progressLabel.setText(f"NCs processed: {processed_ncs}/{nc_total} {self.percentage:.3f}%")
        self.organismProgressLabel.setText(f"Organisms processed: {processed_organism}/{total_organism} {self.organismPercentage:.3f}%")
        self.progressBar.setValue(int(self.percentage))
        self.progressBar.setFormat(f"{self.percentage:.3f}%")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    mainWindow = MainWindow()
    mainWindow.show()
    sys.exit(app.exec_())
