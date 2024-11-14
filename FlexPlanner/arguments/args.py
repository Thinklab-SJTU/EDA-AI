import argparse

def get_args():

    parser = argparse.ArgumentParser(description='Run the 3D floorplanning algorithm')
    parser.add_argument('--lr', type=float, default=1e-4, help='learning rate')
    parser.add_argument('--overlap_ratio', '-o', type=float, default=0.1, help='overlap ratio')
    parser.add_argument('--along_boundary', '-b', type=int, default=1, help='along boundary')
    parser.add_argument('--num_grid_x', type=int, default=64, help='number of grid x')
    parser.add_argument('--num_grid_y', type=int, default=64, help='number of grid y')
    parser.add_argument('--circuit', '-c', type=str, default="n10", help='circuit name')
    parser.add_argument('--result_dir', '-r', type=str, default="result-debug", help='result directory')
    parser.add_argument('--area_util', '-u', type=float, default=1.6, help='area utilization')
    parser.add_argument('--seed', type=int, default=3407, help='random seed')
    parser.add_argument('--num_env', '-e', type=int, default=16, help='num of env')
    parser.add_argument('--num_env_test', type=int, default=4, help='num of envfor test')
    parser.add_argument('--wiremask_bbo', type=int, default=0, help='wiremask bbo')
    parser.add_argument('--device', '-d', type=str, default='cpu', help='device')
    parser.add_argument('--max_epoch', type=int, default=1000, help='max epoch')
    parser.add_argument('--repeat_per_collect', type=int, default=10, help='repeat per collect')
    parser.add_argument('--batch_size', type=int, default=64, help='batch size')
    parser.add_argument('--shared_encoder_final_shape', nargs='+', type=int, default=[4,16,16], help='shared_encoder_final_shape')
    parser.add_argument('--episode_per_collect_per_env', type=int, default=1, help='episode per collect per env')
    parser.add_argument('--save_fig', type=int, default=20, help='save fig gap, <=0 means disable')
    parser.add_argument('--save_checkpoint_interval', type=int, default=50, help='save checkpoint interval')
    parser.add_argument('--train', type=int, default=1, help='train')
    parser.add_argument('--save_batch', type=int, default=0, help='load checkpoint')
    parser.add_argument('--num_preplaced_module', type=int, default=0, help='number of pre-placed module each die.')
    parser.add_argument('--place_order_sorting_method', type=str, default='area', help='place order sorting method')
    parser.add_argument('--add_virtual_block', type=int, default=0, help='the first block to place')
    parser.add_argument('--num_layer', type=int, default=2, help='number of layers/dies')

    # checkpoint
    parser.add_argument('--checkpoint', type=str, default=None, help='checkpoint')
    parser.add_argument('--load_then_collect', type=int, default=1, help='load then collect')
    parser.add_argument('--load_optimizer', type=int, default=0, help='load optimizer, continue on last state, restore')
    parser.add_argument('--statistics', type=str, default=None, help='path to statistics for normalization')
    parser.add_argument('--statistics_method', type=int, default=0, help='statistics method')

    # for block alignment
    parser.add_argument('--enable_alignment', type=int, default=1, help='upper die and bottom die should be aligned')
    parser.add_argument('--num_alignment', type=int, default=10, help='num of alignment pair')
    parser.add_argument('--alignment_rate', type=float, default=0.1, help='c * min(a1, a2) as the threshold.')
    parser.add_argument('--alignment_sort', type=str, default="area", help='alignment sorting method')
    parser.add_argument('--input_alignment_mask', type=int, default=1, help='input alignment mask to AC')
    parser.add_argument('--use_alignment_constraint', type=int, default=1, help='whether use alignment mask to filter out invalid actions')
    parser.add_argument('--reward_weight_alignment', type=float, default=0.05, help='reward weight alignment')


    # for some ablation study
    parser.add_argument('--max_grad_norm', type=float, default=-1.0, help='max grad norm, < 0 means no grad norm')
    parser.add_argument('--norm_wiremask', type=int, default=0, help='norm wiremask')
    parser.add_argument('--place_order_die_by_die', type=int, default=0, help='place order die by die')
    parser.add_argument('--set_vision_to_zero', type=int, default=0, help='set vision to zero')
    parser.add_argument('--set_canvas_to_zero', type=int, default=0, help='set canvas to zero')

    # for all visual input (all masks)
    parser.add_argument('--input_next_block', type=int, default=1, help='input wiremask and position mask of next block', choices=[0, 1])
    parser.add_argument('--input_partner_die', type=int, default=1, help='input partner die information to AC')


    # for shared encoder
    parser.add_argument('--shared_encoder_cls', type=str, default='SharedEncoder', help='shared encoder class')
    parser.add_argument('--graph', type=int, default=0, help='input graph into shared encoder')


    # for (position) actor
    parser.add_argument('--actor_update_shared_encoder', type=int, default=0, help='whether update shared encoder in actor or not')
    parser.add_argument('--deconv_class', type=str, default='InfoGANGenerator', help='deconv class')


    # for critic network
    parser.add_argument('--input_sequence_critic', type=str, default=None, help='input sequence into critic')
    parser.add_argument('--input_die_critic', type=int, default=0, help='input die into critic')
    parser.add_argument('--reduced_dim_critic', type=int, default=0, help='the reduced dim of shared encoder in critic. If <= 0, use the original dim')


    # for block ratio
    parser.add_argument('--enable_ratio', type=int, default=1, help='enable change ratio of movable block')
    parser.add_argument('--ratio_range', nargs='+', type=float, default=[0.75, 1.333], help='ratio range, [min, max]')
    parser.add_argument('--ratio_coef', type=float, default=0.01, help='ratio coef in clip loss')
    parser.add_argument('--ratio_area_in_dim', type=int, default=0, help='ratio decider area dim')
    parser.add_argument('--ratio_share_with_critics', type=int, default=0, help='ratio whether use output from shared encoder or not')


    # for reward function
    parser.add_argument('--reward_func', type=int, default=0, help='reward function')
    parser.add_argument('--use_last_step_reward_to_replace_other_steps', type=int, default=0, help='use last step reward to replace other steps')
    parser.add_argument('--add_last_step_reward_to_other_steps', type=int, default=0, help='add last step reward to other steps')
    parser.add_argument('--reward_weight_hpwl', type=float, default=1.0, help='reward weight hpwl')
    parser.add_argument('--reward_weight_overlap', type=float, default=0.5, help='reward weight overlap')
    parser.add_argument('--reward_weight_final_hpwl', type=float, default=0.0, help='reward weight final hpwl')

    # for original clip loss
    parser.add_argument('--pos_coef', type=float, default=1.0, help='position coef in clip loss')
    parser.add_argument('--ent_coef', type=float, default=0.01, help='entropy coefficient for all entropy loss')


    # for async place
    parser.add_argument('--async_place', type=int, default=1, help='async place')
    parser.add_argument('--async_place_coef', type=float, default=1.0, help='async place coef in clip loss')
    parser.add_argument('--async_place_entropy_coef', type=float, default=1.0, help='async place entropy coef')
    parser.add_argument('--async_place_share_with_critics', type=int, default=0, help='async place whether use output from shared encoder or not')
    parser.add_argument('--async_place_input_sequence', type=str, default=None, help='class of input sequence to layer decider')
    parser.add_argument('--die_embedding', type=int, default=0, help='die embedding for layer decision')
    parser.add_argument('--layer_decider_update_gap', type=int, default=1, help='layer decider update gap')
    parser.add_argument('--learn_async_in_sync', type=int, default=None, help='Deprecated')
    parser.add_argument('--input_layer_sequence', type=int, default=0, help='input layer sequence to layer decider')

    args = parser.parse_args()

    args.wiremask_bbo = True if args.wiremask_bbo > 0 else False
    args.along_boundary = True if args.along_boundary > 0 else False
    args.max_grad_norm = None if args.max_grad_norm < 0 else args.max_grad_norm
    args.actor_update_shared_encoder = True if args.actor_update_shared_encoder > 0 else False
    args.shared_encoder_final_shape = tuple(args.shared_encoder_final_shape)
    args.place_order_die_by_die = True if args.place_order_die_by_die > 0 else False
    args.input_sequence_critic = None if isinstance(args.input_sequence_critic, str) and args.input_sequence_critic.lower() == 'none' else args.input_sequence_critic
    args.save_batch = True if args.save_batch > 0 else False
    args.use_last_step_reward_to_replace_other_steps = True if args.use_last_step_reward_to_replace_other_steps > 0 else False
    args.add_last_step_reward_to_other_steps = True if args.add_last_step_reward_to_other_steps > 0 else False
    args.load_then_collect = True if args.load_then_collect > 0 else False
    args.add_virtual_block = True if args.add_virtual_block > 0 else False

    # override the default value for alignment
    args.enable_alignment = True if args.enable_alignment > 0 else False
    args.num_alignment = args.num_alignment if args.enable_alignment else None
    args.alignment_rate = args.alignment_rate if args.enable_alignment else None
    args.reward_weight_alignment = args.reward_weight_alignment if args.enable_alignment else None
    args.alignment_sort = args.alignment_sort if args.enable_alignment else None
    if args.enable_alignment:
        args.use_alignment_constraint = True if args.use_alignment_constraint > 0 else False
        args.input_alignment_mask = True if args.input_alignment_mask > 0 else False
    else:
        args.use_alignment_constraint = None
        args.input_alignment_mask = None


    # override the default value for change block aspect ratio
    args.enable_ratio = True if args.enable_ratio > 0 and not args.wiremask_bbo else False
    args.ratio_range = args.ratio_range if args.enable_ratio else None
    args.ratio_coef = args.ratio_coef if args.enable_ratio else None
    args.ratio_area_in_dim = args.ratio_area_in_dim if args.enable_ratio else None
    args.ratio_share_with_critics = True if args.ratio_share_with_critics > 0 and args.enable_ratio else False

    # override the default value for async place
    args.async_place = True if args.async_place > 0 and not args.wiremask_bbo else False
    args.async_place_coef = args.async_place_coef if args.async_place else None
    args.async_place_entropy_coef = args.async_place_entropy_coef if args.async_place else None

    args.async_place_share_with_critics = True if args.async_place_share_with_critics > 0 and args.async_place else False
    args.die_embedding = True if args.die_embedding > 0 and args.async_place else False
    args.input_layer_sequence = True if args.input_layer_sequence > 0 and args.async_place else False
    args.async_place_input_sequence = args.async_place_input_sequence if args.async_place else None

    if isinstance(args.async_place_input_sequence, str) and args.async_place_input_sequence.lower() == 'none':
        args.async_place_input_sequence = None
    

    args.train = True if args.train > 0 else False
    args.load_optimizer = True if args.checkpoint is not None and args.load_optimizer > 0 and args.train else False
    args.norm_wiremask = True if args.norm_wiremask > 0 else False
    args.input_partner_die = True if args.input_partner_die > 0 else False

    return args