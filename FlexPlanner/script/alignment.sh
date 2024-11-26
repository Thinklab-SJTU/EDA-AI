export PYTHONUNBUFFERED=1
export CUDA_VISIBLE_DEVICES=0


device=cuda
num_env=8
num_env_test=2
overlap_ratio=0.1
ent_coef=0.01
num_preplaced_module=0
num_layer=2
add_virtual_block=1


circuit=n100
max_epoch=2000
area_util=1.6
num_grid=128
train=1
save_fig=10
lr=1e-4
batch_size=128


# alignment settings
enable_alignment=1
num_alignment=30
alignment_rate=0.1
alignment_sort=area
input_alignment_mask=1
use_alignment_constraint=1
reward_weight_alignment=0.5

shared_encoder_cls=SharedEncoder

graph=2

reward_func=5
reward_weight_hpwl=1.0
reward_weight_overlap=0.5
reward_weight_final_hpwl=1.0


use_last_step_reward_to_replace_other_steps=0
add_last_step_reward_to_other_steps=1


# input_sequence_critic=None
# input_sequence_critic=CNN
# input_sequence_critic=CNN2
input_sequence_critic=Transformer2


input_die_critic=1
input_next_block=1
input_partner_die=1

# learn aspect ratio for blocks
enable_ratio=1
ratio_coef=0.01
ratio_area_in_dim=0
ratio_share_with_critics=0
ratio_min=0.75
ratio_max=1.333

# asynchrnous layer decision
async_place=0
async_place_coef=1.0
async_place_entropy_coef=1.0
async_place_share_with_critics=0
die_embedding=1
# async_place_input_sequence=None
# async_place_input_sequence=CNN
# async_place_input_sequence=CNN2
async_place_input_sequence=Transformer2


place_order_die_by_die=0
norm_wiremask=0
statistics_method=0
load_then_collect=1

# [1] resume training: [checkpoint], [load_optimizer], and [statistics]
# [2] continue training for more epochs: [checkpoint], [load_optimizer], [statistics] and [max_epoch_new]
# [3] transfer learning: [checkpoint]
# [4] test: [checkpoint]
checkpoint=
load_optimizer=
statistics=
max_epoch_new=


if [ -n "${max_epoch_new}" ]; then
    max_epoch=${max_epoch_new}
fi


additional_args=
if [ -n "${checkpoint}" ]; then
    additional_args="${additional_args} --checkpoint ${checkpoint}"
fi
if [ -n "${load_optimizer}" ]; then
    additional_args="${additional_args} --load_optimizer ${load_optimizer}"
fi
if [ -n "${statistics}" ]; then
    additional_args="${additional_args} --statistics ${statistics}"
fi


result_dir=result/${circuit}


# if have defined the variable `checkpoint`, add more comment to the `result_dir`
if [ -n "${checkpoint}" ]; then
    if [ ${train} -eq 1 ]; then
        # resume training or fine-tune
        result_dir=${result_dir}-resume
    else
        # test
        result_dir=${result_dir}-test
    fi
    result_dir=${result_dir}-$(date "+%Y%m%d-%H%M%S")
fi


# if result dir does not exist, create it, else remove it and create a new one
if [ -d ${result_dir} ]; then
    rm -rf ${result_dir}
fi
mkdir -p ${result_dir}


# copy this file to result_dir
cp $0 ${result_dir}/command.sh
touch ${result_dir}/$(hostname)

# copy the checkpoint to result_dir
if [ -n "${checkpoint}" ]; then
    cp ${checkpoint} ${result_dir}/finetune_base_checkpoint.pt
fi

nohup python -u main.py \
    --circuit ${circuit} \
    --async_place ${async_place} \
    --result_dir ${result_dir} \
    --num_env ${num_env} \
    --device ${device} \
    --max_epoch ${max_epoch} \
    --shared_encoder_final_shape 4 16 16 \
    --ratio_range ${ratio_min} ${ratio_max} \
    --ent_coef ${ent_coef} \
    --area_util ${area_util} \
    --overlap_ratio ${overlap_ratio} \
    --enable_ratio ${enable_ratio} \
    --ratio_coef ${ratio_coef} \
    --ratio_area_in_dim ${ratio_area_in_dim} \
    --ratio_share_with_critics ${ratio_share_with_critics} \
    --async_place_coef ${async_place_coef} \
    --async_place_entropy_coef ${async_place_entropy_coef} \
    --async_place_share_with_critics ${async_place_share_with_critics} \
    --enable_alignment ${enable_alignment} \
    --num_alignment ${num_alignment} \
    --alignment_rate ${alignment_rate} \
    --norm_wiremask ${norm_wiremask} \
    --input_partner_die ${input_partner_die} \
    --input_alignment_mask ${input_alignment_mask} \
    --shared_encoder_cls ${shared_encoder_cls} \
    --alignment_sort ${alignment_sort} \
    --die_embedding ${die_embedding} \
    --input_next_block ${input_next_block} \
    --place_order_die_by_die ${place_order_die_by_die} \
    --async_place_input_sequence ${async_place_input_sequence} \
    --input_sequence_critic ${input_sequence_critic} \
    --input_die_critic ${input_die_critic} \
    --reward_func ${reward_func} \
    --graph ${graph} \
    --use_last_step_reward_to_replace_other_steps ${use_last_step_reward_to_replace_other_steps} \
    --add_last_step_reward_to_other_steps ${add_last_step_reward_to_other_steps} \
    --reward_weight_hpwl ${reward_weight_hpwl} \
    --reward_weight_overlap ${reward_weight_overlap} \
    --reward_weight_alignment ${reward_weight_alignment} \
    --reward_weight_final_hpwl ${reward_weight_final_hpwl} \
    --statistics_method ${statistics_method} \
    --load_then_collect ${load_then_collect} \
    --num_grid_x ${num_grid} \
    --num_grid_y ${num_grid} \
    --train ${train} \
    --save_fig ${save_fig} \
    --num_env_test ${num_env_test} \
    --num_layer ${num_layer} \
    --num_preplaced_module ${num_preplaced_module} \
    --lr ${lr} \
    --batch_size ${batch_size} \
    --add_virtual_block ${add_virtual_block} \
    --use_alignment_constraint ${use_alignment_constraint} \
    ${additional_args} \
    > ${result_dir}/log.txt 2>&1 &


