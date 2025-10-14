#!/usr/bin/env python

import sys
from typing import List, Dict
from PyQt5.QtWidgets import QApplication, QWidget, QLabel, QVBoxLayout, QHBoxLayout, QLineEdit, QPushButton, QSizePolicy, QFileDialog, QSpinBox, QComboBox, QGridLayout, QCheckBox, QDoubleSpinBox, QFrame, QRadioButton, QTextEdit, QAbstractButton, QWidget, QMessageBox
from PyQt5.QtCore import Qt, QProcess


#QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)

app: QApplication = QApplication(sys.argv)

class FileEditBox(QFrame):
    def __init__(self, caption: str, filter: str, load: bool):
        super(FileEditBox, self).__init__()

        self.textBox = QLineEdit()
        self.button = QPushButton("...")

        self.caption = caption
        self.filter = filter
        self.load = load
        
        self.textBox.setMinimumWidth(300)
        self.button.clicked.connect(self.onButtonClicked)
        
        layout = QHBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.textBox)
        layout.addWidget(self.button)
        self.setLayout(layout)

    def onButtonClicked(self) -> None:
        if self.load:
            name, type = QFileDialog.getOpenFileName(window, self.caption, self.text(), self.filter, self.filter)
        else:
            name, type = QFileDialog.getSaveFileName(window, self.caption, self.text(), self.filter, self.filter)
        if name != "":
            self.textBox.setText(name)

    def text(self) -> str:
        return self.textBox.text()


def createAminoCheckBox(name: str, letters: str, commandLineOption: str) -> QCheckBox:
    check = QCheckBox()
    check.letters = letters
    check.setText(name)
    check.commandLineOption = commandLineOption
    return check

aminoCheckBoxes : List[QCheckBox] = [
    createAminoCheckBox("Aliphatic", "IVL", "-aliphatic"),
    createAminoCheckBox("Sulphur", "MC", "-sulphur"),
    createAminoCheckBox("Tiny", "AGCS", "-tiny"),
    createAminoCheckBox("Aromatic", "FYWH", "-aromatic"),
    createAminoCheckBox("Hydrophobic", "ACTKHWYFMILV", "-hydrophobic"),
    createAminoCheckBox("Charged", "DEHKR", "-charged"),
    createAminoCheckBox("Positive", "HKR", "-positiv"),
    createAminoCheckBox("Polar", "NQSTCDEHKRYW", "-polar"),
    createAminoCheckBox("Acidic", "NQ", "-acidic"),
    createAminoCheckBox("Small", "VPAGCTSDN", "-small"),
    createAminoCheckBox("Hydroxylic", "ST", "-hydroxylic")
]

def onProteinWildcardClicked(checked):
    letters = set()
    for cb in aminoCheckBoxes:
        if cb.isChecked():
            for letter in cb.letters:
                letters.add(letter)
    for cb in aminoCheckBoxes:
        if not cb.isChecked():
            cb.setEnabled(letters.isdisjoint(cb.letters))

def bindToggleButtonToEnablements(button: QAbstractButton, control: QWidget):
    listener = lambda checked: control.setEnabled(checked)
    button.toggled.connect(listener)
    control.setEnabled(button.isChecked())

def withMargin(widget: QWidget, left: int, top: int, right: int, bottom: int) -> QHBoxLayout:
    box = QHBoxLayout()
    box.addWidget(widget)
    box.setContentsMargins(left, top, right, bottom)
    return box

dnaRnaRepeatTypes: List[str] = ["direct non-complementary repeats", "inverse non-complementary repeats", "direct complementary repeats", "inverse complementary repeats"]
proteinRepeatTypes: List[str] = ["direct non-complementary repeats", "inverse non-complementary repeats"]

sequenceTypeTextCommandLineSwitch: Dict[str, str] = {
    "DNA": "-dna",
    "RNA": "-rna",
    "Protein": "-protein"
}

repeatTypesTextToCommandLineSwitch: Dict[str, str] = {
    "direct non-complementary repeats": "-dn",
    "inverse non-complementary repeats": "-in",
    "direct complementary repeats": "-dc",
     "inverse complementary repeats": "ic"
}

gridRow: int = 0
grid: QGridLayout = QGridLayout()

documentation: QLabel = QLabel()
documentation.setText('<a href="RepeatsHelp.html">Documentation</a>')
documentation.setOpenExternalLinks(True)
grid.addWidget(documentation, gridRow, 0)
gridRow += 1

grid.addWidget(QLabel("FASTA file:"), gridRow, 0)
fileNameBox: FileEditBox = FileEditBox("Open FASTA file", "*.fasta", True)
grid.addWidget(fileNameBox, gridRow, 1, 1, 2)
gridRow += 1

grid.addWidget(QLabel("Sequence type: "), gridRow, 0)
sequenceTypeCombo: QComboBox = QComboBox()
sequenceTypeCombo.addItems(["DNA", "RNA", "Protein"])
grid.addWidget(sequenceTypeCombo, gridRow, 1, 1, 2, Qt.AlignLeft)
gridRow += 1

grid.addWidget(QLabel("Minimal repeat length: "), gridRow, 0)
minRepeatLengthSpinBox: QSpinBox = QSpinBox()
minRepeatLengthSpinBox.setMinimum(1)
minRepeatLengthSpinBox.setValue(8)
minRepeatLengthSpinBox.setMinimumWidth(60)
grid.addWidget(minRepeatLengthSpinBox, gridRow, 1, 1, 2, Qt.AlignLeft)
gridRow += 1

grid.addWidget(QLabel("Repeat type: "), gridRow, 0)
repeatTypeCombo: QComboBox = QComboBox()
repeatTypeCombo.addItems(dnaRnaRepeatTypes)
grid.addWidget(repeatTypeCombo, gridRow, 1, 1, 2, Qt.AlignLeft)
gridRow += 1

maxGapCheckBox: QCheckBox = QCheckBox()
maxGapCheckBox.setText("Max gap:")
grid.addWidget(maxGapCheckBox, gridRow, 0)
maxGapSpinBox: QDoubleSpinBox = QSpinBox()
maxGapSpinBox.setMinimum(1)
maxGapSpinBox.setMinimumWidth(60)
grid.addWidget(maxGapSpinBox, gridRow, 1, 1, 2, Qt.AlignLeft)
bindToggleButtonToEnablements(maxGapCheckBox, maxGapSpinBox)
gridRow += 1

ambiguousLettersCheckBox: QCheckBox = QCheckBox()
ambiguousLettersCheckBox.setText("Ambiguous letters")
grid.addWidget(ambiguousLettersCheckBox, gridRow, 0, 1, 2)
gridRow += 1

motifCheckBox: QCheckBox = QCheckBox()
motifCheckBox.setText("Motif")
grid.addWidget(motifCheckBox, gridRow, 0)
motifText: QLineEdit() = QLineEdit()
grid.addWidget(motifText, gridRow, 1)
bindToggleButtonToEnablements(motifCheckBox, motifText)
motifMismatchCheckBox: QCheckBox = QCheckBox()
motifMismatchCheckBox.setText("Mismatches")
grid.addWidget(motifMismatchCheckBox, gridRow, 2)
bindToggleButtonToEnablements(motifCheckBox, motifMismatchCheckBox)
gridRow += 1

useProbabilityCheckBox: QCheckBox = QCheckBox()
useProbabilityCheckBox.setChecked(True)
useProbabilityCheckBox.setText("P value:")
grid.addWidget(useProbabilityCheckBox, gridRow, 0)
pValueSpinBox: QDoubleSpinBox = QDoubleSpinBox()
pValueSpinBox.setDecimals(2)
pValueSpinBox.setMinimum(0.01)
pValueSpinBox.setMaximum(0.99)
pValueSpinBox.setSingleStep(0.01)
pValueSpinBox.setValue(0.95)
pValueSpinBox.setMinimumWidth(60)
grid.addWidget(pValueSpinBox, gridRow, 1, 1, 2, Qt.AlignLeft)
bindToggleButtonToEnablements(useProbabilityCheckBox, pValueSpinBox)
gridRow += 1

useCustomLetterMapping: QCheckBox = QCheckBox("Use a custom complement mapping file:")
grid.addWidget(useCustomLetterMapping, gridRow, 0)
customLetterMappingFile: FileEditBox = FileEditBox("Open complement mapping file", "*.complements", True)
grid.addWidget(customLetterMappingFile, gridRow, 1, 1, 2)
bindToggleButtonToEnablements(useCustomLetterMapping, customLetterMappingFile)
gridRow += 1

proteinWildcardsLabel: QLabel = QLabel("Protein wildcards:")
grid.addWidget(proteinWildcardsLabel, gridRow, 0, Qt.AlignTop)
wildcardGrid: QGridLayout = QGridLayout()
wildcardGrid.setContentsMargins(0, 0, 0, 0)
for aminoPos, check in enumerate(aminoCheckBoxes):
    check.clicked.connect(onProteinWildcardClicked)
    aminosColumns = 4
    wildcardGrid.addWidget(check, aminoPos / aminosColumns, aminoPos % aminosColumns)
wildcardFrame: QFrame = QFrame()
wildcardFrame.setContentsMargins(0, 0, 0, 0)
wildcardFrame.setLayout(wildcardGrid)
wildcardFrame.setEnabled(False)
grid.addWidget(wildcardFrame, gridRow, 1, 1, 2)
gridRow += 1

multiSequence: QCheckBox = QCheckBox()
multiSequence.setText("Treat all sequences in the FASTA file as a single sequence")
grid.addLayout(withMargin(multiSequence, 0, 10, 0, 0), gridRow, 0, 1, 2)
gridRow += 1

ignoreN: QCheckBox = QCheckBox();
ignoreN.setText("Exclude letter N segments length:")
grid.addWidget(ignoreN, gridRow, 0)
ignoreNLen: QSpinBox = QSpinBox()
ignoreNLen.setMinimum(1)
ignoreNLen.setMinimumWidth(60)
grid.addWidget(ignoreNLen, gridRow, 1, 1, 2, Qt.AlignLeft)
bindToggleButtonToEnablements(ignoreN, ignoreNLen)
gridRow += 1

grid.addLayout(withMargin (QLabel("Output:"), 0, 20, 0, 0), gridRow, 0)
gridRow += 1

fileOutputRadio: QRadioButton = QRadioButton("Output file name:")
grid.addWidget(fileOutputRadio, gridRow, 0)
fileOutputDestination: FileEditBox = FileEditBox("Output file", "*.txt", False)
fileOutputDestination.setEnabled(False)
grid.addWidget(fileOutputDestination, gridRow, 1, 1, 2)
bindToggleButtonToEnablements(fileOutputRadio, fileOutputDestination)
gridRow += 1

databaseOutputRadio: QRadioButton = QRadioButton("ODBC connection string:")
databaseConnectionString: QLineEdit = QLineEdit()
grid.addWidget(databaseOutputRadio, gridRow, 0)
grid.addWidget(databaseConnectionString, gridRow, 1, 1, 2)
bindToggleButtonToEnablements(databaseOutputRadio, databaseConnectionString)
gridRow += 1

inputsFrame: QFrame = QFrame();
inputsFrame.setContentsMargins(0, 0, 0, 0)
inputsFrame.setLayout(grid)

vbox: QVBoxLayout = QVBoxLayout()
vbox.addWidget(inputsFrame)

grid = QGridLayout()
gridRow = 0

runButton: QPushButton = QPushButton("Start")
grid.addWidget(runButton, gridRow, 0, 1, 2, Qt.AlignRight)
gridRow += 1

consoleOutput: QTextEdit = QTextEdit()
consoleOutput.setReadOnly(True)
grid.addWidget(consoleOutput, gridRow, 0, 1, 2)
gridRow += 1

vbox.addLayout(grid)

window: QWidget = QWidget()
window.setLayout(vbox)
window.setWindowTitle("RepeatPlus")  # Here we set the title for our window.


def sequenceTypeComboChanged():
    wildcardFrame.setEnabled(sequenceTypeCombo.currentIndex() == 2)
    previousIndex = repeatTypeCombo.currentIndex()
    repeatTypeCombo.clear()
    if sequenceTypeCombo.currentIndex() == 2:
        repeatTypeCombo.addItems(proteinRepeatTypes)
        if (previousIndex < len(proteinRepeatTypes)):
            repeatTypeCombo.setCurrentIndex(previousIndex)
    else:
        repeatTypeCombo.addItems(dnaRnaRepeatTypes)
        repeatTypeCombo.setCurrentIndex(previousIndex)

sequenceTypeCombo.currentIndexChanged.connect(sequenceTypeComboChanged)

fileOutputRadio.setChecked(True)

process: QProcess = QProcess()

def onProcessHasOutput() -> None:
    consoleOutput.setTextColor(Qt.black)
    consoleOutput.append(process.readAllStandardOutput().data().decode('utf8'))

def onProcessHasError() -> None:
    consoleOutput.setTextColor(Qt.red)
    consoleOutput.append(process.readAllStandardError().data().decode('utf8'))

def showError(text: str) -> None:
    QMessageBox.critical(window, "Bad input", text)

def onRunFinished(exitCode: int, exitStatus) -> None:
    runButton.setText("Start")
    inputsFrame.setEnabled(True)

def onAmbiguousOrMotifToggled(checked: bool):
    useProbabilityCheckBox.setEnabled(not ambiguousLettersCheckBox.isChecked() and (not motifCheckBox.isChecked() or not motifMismatchCheckBox.isChecked()))
    pValueSpinBox.setEnabled(useProbabilityCheckBox.isEnabled() and useProbabilityCheckBox.isChecked())

    useCustomLetterMapping.setEnabled(not checked)
    customLetterMappingFile.setEnabled(not checked and useCustomLetterMapping.isChecked())
    proteinWildcardsLabel.setEnabled(not checked)
    wildcardFrame.setEnabled(not checked and sequenceTypeCombo.currentIndex() == 2)
    multiSequence.setEnabled(not checked)
    ignoreN.setEnabled(not checked)
    ignoreNLen.setEnabled(not checked and ignoreN.isChecked())
    if checked:
        fileOutputRadio.setChecked(True)
    databaseOutputRadio.setEnabled(not checked)


def onAmbiguousLettersCheckBoxToggled(checked: bool):
    onAmbiguousOrMotifToggled(checked)
    motifCheckBox.setEnabled(not checked)
    motifText.setEnabled(not checked)

def onMotifCheckBoxToggled(checked: bool):
    onAmbiguousOrMotifToggled(checked)
    ambiguousLettersCheckBox.setEnabled(not checked)

def onMotifMismatchCheckBox(checked: bool):
    useProbabilityCheckBox.setEnabled(not ambiguousLettersCheckBox.isChecked() and (not motifCheckBox.isChecked() or not motifMismatchCheckBox.isChecked()))
    pValueSpinBox.setEnabled(useProbabilityCheckBox.isEnabled() and useProbabilityCheckBox.isChecked())

ambiguousLettersCheckBox.toggled.connect(onAmbiguousLettersCheckBoxToggled)
motifCheckBox.toggled.connect(onMotifCheckBoxToggled)
motifMismatchCheckBox.toggled.connect(onMotifMismatchCheckBox)

def onRunButtonClicked() -> None:
    if inputsFrame.isEnabled():
        fileName: str = fileNameBox.text()
        if fileName == "":
            showError("Please select the input FASTA file.")
            return

        args: List[str] = [
            fileName,
            str(minRepeatLengthSpinBox.value()),
            sequenceTypeTextCommandLineSwitch [sequenceTypeCombo.currentText()],
            repeatTypesTextToCommandLineSwitch [repeatTypeCombo.currentText()]
        ]

        if maxGapCheckBox.isChecked():
            args.append("-maxgap")
            args.append(str(maxGapSpinBox.value()))

        if not (ambiguousLettersCheckBox.isChecked() or (motifCheckBox.isChecked() and motifMismatchCheckBox.isChecked())):
            if useProbabilityCheckBox.isChecked():
                args.append("-pvalue")
                args.append(str(pValueSpinBox.value()))
            else:
                args.append("-noprobability")

        if ambiguousLettersCheckBox.isChecked() or motifCheckBox.isChecked():
            if useCustomLetterMapping.isChecked():
                args.append("-complement")
                complementFileName = customLetterMappingFile.text()
                if (complementFileName == ""):
                    showError("Complement file name is empty.")
                    return;
                args.append(complementFileName)

            if sequenceTypeCombo.currentText() == "Protein":
                for cb in aminoCheckBoxes:
                    if cb.isChecked():
                        args.append(cb.commandLineOption)

            if multiSequence.isChecked():
                args.append("-msr")

            if ignoreN.isChecked():
                args.append("-ex")
                args.append(str(ignoreNLen.value()))

        if motifCheckBox.isChecked():
            if motifText.text() == "":
                showError("Motif is missing")
                return
            args.append("-motif")
            args.append(motifText.text())

        if fileOutputRadio.isChecked():
            args.append("-output")
            if fileOutputDestination.text() == "":
                showError("Output file name is empty.")
                return
            args.append(fileOutputDestination.text())

        if databaseOutputRadio.isChecked():
            args.append("-database")
            if databaseConnectionString.text() == "":
                showError("Database connection string is empty")
                return
            args.append(databaseConnectionString.text())

        consoleOutput.clear()
        runButton.setText("Stop")
        inputsFrame.setEnabled(False)

        processName: str = "StatRepeats"
        if ambiguousLettersCheckBox.isChecked() or (motifCheckBox.isChecked() and motifMismatchCheckBox.isChecked()):
            processName = "Repeats"

        global process
        process = QProcess()
        process.readyReadStandardOutput.connect(onProcessHasOutput)
        process.readyReadStandardError.connect(onProcessHasError)
        process.finished.connect(onRunFinished)
        process.start(processName, args)
    else:
        process.kill()
        runButton.setText("Start")
        inputsFrame.setEnabled(True)


runButton.clicked.connect(onRunButtonClicked)

window.show()  # The show() method displays the widget on the screen.

sys.exit(app.exec_())  # Finally, we enter the mainloop of the application.
