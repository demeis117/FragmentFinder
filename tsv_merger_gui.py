# tsv_merger_gui_check.py
# GUI for Merge_stats_tolerant_check.py
# - Top: "Select Directory with .tsv Files" + small label showing the chosen folder
# - Middle: "Select FASTA Database File" + small label showing the chosen file
# - Bottom: "Run Merge" + small label showing the output (_checked.tsv) after completion

import os, sys, time, traceback
from pathlib import Path

# Use PySide6; if unavailable, fall back to PyQt5
from PySide6 import QtWidgets, QtCore
Signal = QtCore.Signal
Slot = QtCore.Slot

# --- import the merge/check module (must be in the same directory) ---
import Merge_stats_tolerant_check as M  # <-- your uploaded script

def run_merge_pipeline(tsv_dir: str, fasta_path: str) -> str:
    """
    Recreates M.main() but parameterized, and returns the final CHECKED TSV path.
    This mirrors the sequence in Merge_stats_tolerant_check.py.
    """
    tsv_dir = M.sanitize_dir(tsv_dir)
    fasta_path = M.sanitize_dir(fasta_path)

    file_list = M.get_files(tsv_dir)
    if not file_list:
        raise FileNotFoundError(f"No eligible .tsv files found in: {tsv_dir}")

    out_file = M.get_new_out_file_name(file_list)   # e.g. "(FFdb_Final_Mar2025)"
    file_count = M.get_list_count(file_list)

    # 1) merge & intermediate grouping
    M.merge_files(tsv_dir)
    M.init_group(tsv_dir)

    # 2) load DB strings and build final grouped table
    hmap = M.get_orig_strings(fasta_path)
    M.get_final_group(tsv_dir, hmap)
    M.export_final_tsv(tsv_dir, file_count, out_file)

    # 3) cross-check origins; this produces *_checked.tsv
    merged_path = f"{tsv_dir}/Joined_Results/{out_file}.tsv"
    checked_out = os.path.splitext(merged_path)[0] + "_checked.tsv"
    M.check_origins_main(merged_path, fasta_path, checked_out)

    # 4) best-effort cleanup of intermediates (same as script)
    for temp in ("files_and_reads.tsv", "grouped_by_peak.tsv", "t1.tsv"):
        p = Path(tsv_dir) / temp
        if p.exists():
            try: p.unlink()
            except OSError: pass

    return checked_out

# ---------------- Worker (keeps UI responsive) ----------------
class Worker(QtCore.QObject):
    finished = Signal(str, float)  # (output_path, seconds)
    error = Signal(str)

    @Slot(str, str)
    def do_work(self, dir_path, fasta_path):
        t0 = time.time()
        try:
            outp = run_merge_pipeline(dir_path, fasta_path)
            self.finished.emit(outp, time.time() - t0)
        except Exception as e:
            self.error.emit(f"{e}\n\n{traceback.format_exc()}")

# ---------------- Main Window ----------------
class Main(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("TSV Merger")

        # Buttons
        self.btn_dir = QtWidgets.QPushButton("Select Directory with .tsv Files")
        self.btn_fa  = QtWidgets.QPushButton("Select FASTA Database File")
        self.btn_run = QtWidgets.QPushButton("Run Merge")

        # Small path labels under each button
        self.lbl_dir = QtWidgets.QLabel(" ")
        self.lbl_fa  = QtWidgets.QLabel(" ")
        self.lbl_out = QtWidgets.QLabel(" ")

        for lbl in (self.lbl_dir, self.lbl_fa, self.lbl_out):
            lbl.setWordWrap(True)
            lbl.setStyleSheet("color:#444; font-size:11px;")

        # Layout to match your screenshot: button then its label
        v = QtWidgets.QVBoxLayout(self)
        v.addWidget(self.btn_dir); v.addWidget(self.lbl_dir); v.addSpacing(10)
        v.addWidget(self.btn_fa);  v.addWidget(self.lbl_fa);  v.addSpacing(14)
        v.addWidget(self.btn_run); v.addWidget(self.lbl_out); v.addStretch(1)

        # Connections
        self.btn_dir.clicked.connect(self.pick_dir)
        self.btn_fa.clicked.connect(self.pick_fasta)
        self.btn_run.clicked.connect(self.run_merge)

        # Thread/worker setup
        self.thread = QtCore.QThread(self)
        self.worker = Worker()
        self.worker.moveToThread(self.thread)
        self.worker.finished.connect(self.on_finished)
        self.worker.error.connect(self.on_error)
        self.thread.start()

        # State
        self.tsv_dir = ""
        self.fasta   = ""

        # Make buttons wide like the mockup
        for b in (self.btn_dir, self.btn_fa, self.btn_run):
            b.setMinimumWidth(480)

    def pick_dir(self):
        d = QtWidgets.QFileDialog.getExistingDirectory(self, "Select directory containing .tsv files")
        if d:
            self.tsv_dir = d
            self.lbl_dir.setText(d)

    def pick_fasta(self):
        fn, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Select FASTA database file",
            filter="FASTA (*.fa *.fasta *.fna *.txt);;All Files (*)"
        )
        if fn:
            self.fasta = fn
            self.lbl_fa.setText(fn)

    def run_merge(self):
        if not self.tsv_dir:
            QtWidgets.QMessageBox.warning(self, "Missing", "Please select the TSV directory first.")
            return
        if not self.fasta:
            QtWidgets.QMessageBox.warning(self, "Missing", "Please select the FASTA database file.")
            return
        self.btn_run.setEnabled(False)
        self.lbl_out.setText("Running…")
        QtCore.QMetaObject.invokeMethod(
            self.worker, "do_work", QtCore.Qt.QueuedConnection,
            QtCore.Q_ARG(str, self.tsv_dir),
            QtCore.Q_ARG(str, self.fasta)
        )

    def on_finished(self, out_path, seconds):
        self.btn_run.setEnabled(True)
        self.lbl_out.setText(f"Output: {out_path}\nElapsed: {seconds:.1f}s")
        QtWidgets.QMessageBox.information(self, "Done", f"Merged (checked) TSV written to:\n{out_path}")

    def on_error(self, msg):
        self.btn_run.setEnabled(True)
        self.lbl_out.setText("Error. See details.")
        QtWidgets.QMessageBox.critical(self, "Error", msg)

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    w = Main()
    w.resize(560, 200)
    w.show()
    sys.exit(app.exec())
