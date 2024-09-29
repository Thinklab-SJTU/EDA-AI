export CUDA_VISIBLE_DEVICES=0
seed=3407
graph_autoencoder_checkpoint="path/to/pre_trained/graph/encoder"
python main.py \
    --root_dir "result-PreRoutGNN" \
    --model_type "PreRoutGNN" \
    --normalize_pin_location \
    --use_graph_autoencoder \
    --graph_autoencoder_checkpoint $graph_autoencoder_checkpoint \
    --seed $seed \
    --device "cuda" \
    --move_to_cuda_in_advance 0 \
    --graph_autoencoder_lr 0.0005 \
    --ema