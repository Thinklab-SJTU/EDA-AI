import os
from options.test_options import TestOptions
from data import create_dataset
from models import create_model
from utils import visualizer

if __name__ == '__main__':
    opt = TestOptions().parse()  # get test options
    # hard-code some parameters for test
    opt.num_threads = 0
    opt.batch_size = 1
    opt.serial_batches = True  # disable data shuffling; comment this line if results on randomly chosen
    # images are needed.
    opt.no_flip = True  # no flip; comment this line if results on flipped images are needed.
    dataset = create_dataset(opt)  # create a dataset given opt.dataset_mode and other options
    dataset_size = len(dataset)
    model = create_model(opt)  # create a model given opt.model and other options
    model.setup(opt)  # regular setup: load and print networks; create schedulers

    # test with eval mode. This only affects layers like batchnorm and dropout.
    if opt.eval:
        model.eval()
    for i, data in enumerate(dataset):
        if i >= opt.num_test:  # only apply our model to opt.num_test images.
            break
        model.set_input(data)  # unpack data from data loader
        model.test()  # run inference
        visuals = model.get_current_visuals()  # get image results
        img_path = model.get_image_paths()  # get image paths
        if i % 250 == 0:  # save images to the disk
            visualizer.display_test_results(os.path.join(opt.results_dir, opt.name + '_' + opt.epoch), visuals, i)

    print('Correctness: %f, TWL Ratio: %f' % (model.correctness / dataset_size, model.aligned_fake_twl / (model.aligned_real_twl + 1e-5)))
