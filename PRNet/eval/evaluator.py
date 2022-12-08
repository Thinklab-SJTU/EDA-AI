import matplotlib.pyplot as plt
import os
import numpy as np

# matplotlib.use('Gdk')


class Evaluator:
    def __init__(self, opt):
        self.path = os.path.join(opt.checkpoints_dir, opt.name)
        self.correctness = []
        self.wl_ratio = []
        self.congestion = []
        self.GAN_loss = []
        self.D_real = []
        self.D_fake = []

    def update(self, correctness, wl_ratio, D_real, D_fake, GAN_loss):
        self.correctness.append(correctness)
        self.wl_ratio.append(wl_ratio)
        self.GAN_loss.append(GAN_loss)
        self.D_real.append(D_real)
        self.D_fake.append(D_fake)
        # self.congestion.append(congestion)

    def plot(self):
        epoch_list = list(range(len(self.correctness)))
        plt.plot(epoch_list, self.correctness)
        plt.xlim(0, len(epoch_list))
        plt.xlabel("Epoch")
        plt.ylabel("Correctness")
        plt.title("Training Correctness")
        plt.savefig(self.path + "/Correctness.jpg")
        plt.clf()

        plt.plot(epoch_list, self.wl_ratio)
        plt.xlim(0, len(epoch_list))
        plt.ylim(0.99, 1.2)
        plt.xlabel("Epoch")
        plt.ylabel("WL Ratio")
        plt.title("Training WireLength Ratio")
        plt.savefig(self.path + "/Wirelength_ratio.jpg")
        plt.clf()

        plt.plot(epoch_list, self.GAN_loss, color='orange', label='GAN loss')
        plt.plot(epoch_list, self.D_real, color='blue', label='D real')
        plt.plot(epoch_list, self.D_fake, color='skyblue', label='D fake')
        plt.xlim(0, len(epoch_list))
        # plt.ylim(-0.1, 1.0)
        plt.title("Training loss")
        plt.xlabel("Epoch")
        plt.ylabel("Loss")
        plt.savefig(self.path + "/Loss.jpg")

        # plt.plot(epoch_list, self.congestion)
        # plt.xlim(0, len(epoch_list))
        # plt.ylim(0, np.median(np.array(self.congestion)) * 2)
        # plt.xlabel("Epoch")
        # plt.ylabel("Congestion")
        # plt.title("Training Congestion")
        # plt.savefig(self.path + "/congestion.jpg")
