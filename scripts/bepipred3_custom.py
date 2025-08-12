### IMPORTS AND STATIC PATHS ###
import os
import sys
import argparse
import zipfile
import shutil
from pathlib import Path

bp3_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src", "BepiPred3.0-Predictor"))
sys.path.insert(0, bp3_path)

from bp3 import bepipred3

WORK_DIR = Path(bp3_path).resolve()

### COMMAND LINE ARGUMENTS ###
parser = argparse.ArgumentParser("Make B-cell epitope predictions from fasta file.")
parser.add_argument("-i", required=True, action="store", dest="fasta_file", type=Path, help="Fasta file contianing antigens")
parser.add_argument("-o", required=True, action="store", dest="out_dir", type=Path, help="Output directory to store B-cell epitope predictions.")
parser.add_argument("-pred", action="store", choices=["mjv_pred", "vt_pred"], required=True, dest="pred", help="Majorty vote ensemble prediction or variable threshold predicition on average ensemble posistive probabilities. ")
parser.add_argument("-add_seq_len", action="store_true", dest="add_seq_len", help="Add sequence lengths to esm-encodings. Default is false. On the web server this option is set to true.")
parser.add_argument("-esm_dir", action="store", default= WORK_DIR / "esm_encodings", dest="esm_dir", type=Path, help="Directory to save esm encodings to. Default is current working directory.")
parser.add_argument("-t", action="store", default=0.1512, type=float, dest="var_threshold", help="Threshold to use, when making predictions on average ensemble positive probability outputs. Default is 0.15.")
parser.add_argument("-top", action="store", default=0.2, type=float, dest="top_cands", help="Top percentage of epitope residues Default is top 20 pct.")
parser.add_argument("-rolling_window_size", default=9, type=int, dest="rolling_window_size", help="Window size to use for rolling average on B-cell epitope probability scores. Default is 9.")
parser.add_argument("-plot_linear_epitope_scores", action="store_true", dest="plot_linear_epitope_scores", help="Use linear B-cell epitope probability scores for plot. Default is false.")
parser.add_argument("-z", action="store_true", dest="zip_results", help="Specify option to create zip the bepipred-3.0 results (except the interactive .html figure). Default is false.")
parser.add_argument("-interactive_plot", action="store_true", dest="interactive_plot",
                    help="Generate interactive HTML plot in addition to other outputs.")
parser.add_argument("-top_pct_files", action="store_true", dest="top_pct_files",
                    help="Generate top percentage epitope residue files. Default is false.")

args = parser.parse_args()
fasta_file = args.fasta_file
out_dir = args.out_dir
var_threshold = args.var_threshold
pred = args.pred
add_seq_len = args.add_seq_len
esm_dir = args.esm_dir
top_cands = args.top_cands
rolling_window_size = args.rolling_window_size
plot_linear_epitope_scores = args.plot_linear_epitope_scores
zip_results = args.zip_results
interactive_plot = args.interactive_plot
top_pct_files = args.top_pct_files

### FUCNCTIONS ###
def zip_function(result_files, outfile):
    zipf = zipfile.ZipFile(outfile, 'w')
    for result_file in result_files: zipf.write(result_file, arcname=result_file.name)
    zipf.close()

### MAIN ###

## Load antigen input and create ESM-2 encodings ## 

#if you have the esm2 model stored locally, you can this command. To work you need both esm2_t33_650M_UR50D.pt and the esm2_t33_650M_UR50D-contact-regression.pt stored in same directory.
#MyAntigens = bepipred3.Antigens(fasta_file, esm_dir, add_seq_len=add_seq_len, run_esm_model_local=str(WORK_DIR / "models" / "esm2_t33_650M_UR50D.pt") )

MyAntigens = bepipred3.Antigens(fasta_file, esm_dir, add_seq_len=add_seq_len)
MyBP3EnsemblePredict = bepipred3.BP3EnsemblePredict(MyAntigens, rolling_window_size=rolling_window_size, top_pred_pct = top_cands)
MyBP3EnsemblePredict.run_bp3_ensemble()

if top_pct_files:
    MyBP3EnsemblePredict.create_toppct_files(out_dir)

## Create predictions .csv  
MyBP3EnsemblePredict.create_csvfile(out_dir)
original_csv_path = out_dir / "raw_output.csv" #by default the output file will be names raw_output.csv 

## Construct new filename based on the fasta filename
fasta_stem = fasta_file.stem  # gets filename without extension
new_csv_name = f"bp3_{fasta_stem}.csv"
new_csv_path = out_dir / new_csv_name

if original_csv_path.exists():
    shutil.copyfile(original_csv_path, new_csv_path)
    print(f"Created {new_csv_path}")
else:
    print(f"Warning: Expected CSV file {original_csv_path} not found.")

## B-cell epitope predictions ##
if pred == "mjv_pred":
    MyBP3EnsemblePredict.bp3_pred_majority_vote(out_dir)
elif pred == "vt_pred":
    MyBP3EnsemblePredict.bp3_pred_variable_threshold(out_dir, var_threshold=var_threshold)

#generate plots (generating graphs for a maximum of 40 proteins)
if args.interactive_plot:
    MyBP3EnsemblePredict.bp3_generate_plots(
        args.out_dir,
        num_interactive_figs=50,
        use_rolling_mean=args.plot_linear_epitope_scores
    )

#zip results
if zip_results:
	print("Zipping results")
	result_files = [f for f in out_dir.glob("*") if f.suffix != ".html"]
	zip_function(result_files, out_dir / "bepipred3_results.zip")
	for result_file in result_files: result_file.unlink()