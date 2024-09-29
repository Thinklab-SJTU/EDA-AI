export CUDA_VISIBLE_DEVICES=0

python main_graph_autoencoder.py \
    --circuits_json_path "circuits_json/circuits-original.json" \
    --data_split_json_path "circuits_split_json/data_split-original.json" \
    --root_dir "result-encoder" \
    --num_level_freq_compents 4 \
    --num_epochs 10000 \
    --device "cuda" \
    --normalize_pin_location