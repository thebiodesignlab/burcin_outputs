#!/bin/bash
# To run ProteinMPNN use this script.
#SBATCH -p gpu
#SBATCH --job-name=1CEZ
#SBATCH --mem=32g
#SBATCH --gres=gpu:1
#SBATCH -c 3
#SBATCH --output=1CEZ.out
#SBATCH --time=4-0:0:0
#SBATCH --mail-user=brcnacar@gmail.com
#SBATCH --mail-type=BEGIN,END

ml Anaconda3/2023.03
conda init bash
source ~/.bashrc
conda activate proteinmpnn

folder_with_pdbs="tempdb"

output_dir="outputs"

if [ ! -d $output_dir ]
then
    mkdir -p $output_dir
fi


path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
path_for_tied_positions=$output_dir"/tied_pdbs.jsonl"
path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"
path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"
chains_to_design="A"
fixed_positions="50 53 54 57 58 60 61 63 68 71 93 94 95 96 97 98 99 100 101 102 104 129 131 133 135 136 139 143 158 160 161 163 164 165 166 167 168 169 171 172 173 174 175 176 177 178 179 182 200 201 206 210 211 215 231 234 235 236 237 238 240 241 242 263 264 265 267 294 295 297 298 299 300 301 304 305 357 377 378 379 380 381 382 385 386 388 389 390 391 392 393 394 396 397 400 419 421 422 423 425 427 429 431 435 436 437 438 441 471 472 527 529 537 538 539 540 541 542 569 571 600 601 602 608 610 627 631 632 635 636 638 639 640 641 642 643 644 645 646 647 648 649 651 652 656 659 660 669 670 671 679 700 704 713 714 722 736 737 738 739 740 741 742 744 745 746 747 748 749 750 753 754 755 756 757 758 759 760 761 762 764 768 771 772 773 775 776 777 780 781 784 787 810 811 812 813 860"

python $pmpnn/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains

python $pmpnn/helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"
python $pmpnn/helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$fixed_positions"

python $pmpnn/helper_scripts/make_pos_neg_tied_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_tied_positions --homooligomer 1 --pos_neg_chain_list="A" --pos_neg_chain_betas "1.0"

python $pmpnn/protein_mpnn_run.py --ca_only \
        --jsonl_path $path_for_parsed_chains \
        --chain_id_jsonl $path_for_assigned_chains \
	--fixed_positions_jsonl $path_for_fixed_positions \
	--tied_positions_jsonl $path_for_tied_positions \
        --out_folder $output_dir \
        --num_seq_per_target 1000 \
        --sampling_temp "0.3" \
	--save_probs 1 \
        --seed 37 --save_score 1 \
        --batch_size 1 \
