export CUDA_VISIBLE_DEVICES=0

pretrained_checkpoint="checkpoint/PreRoutGNN.pt"
graph_autoencoder_checkpoint=$pretrained_checkpoint

python main.py \
    --root_dir "result-inference" \
    --model_type "PreRoutGNN" \
    --num_level_freq_compents 4 \
    --max_level 200 \
    --num_pin_location_freq_compents 0 \
    --normalize_pin_location \
    --use_graph_autoencoder \
    --graph_autoencoder_checkpoint $graph_autoencoder_checkpoint \
    --device "cuda" \
    --test \
    --save_prediction \
    --pretrained_checkpoint $pretrained_checkpoint