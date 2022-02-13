from rtolib.core import NoiseGenerator
import copy

def for_old_version():
    noise_filename="noise/noise1.txt"

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)
    keys = [
                'XFr_A',
                'XFr_B',
                'XFr_E',
                'XFr_P',
                'XFr_G',
        ]

    with open("noise/noise_transformed.txt",'w') as fp:
        for iter in range(30):
            for trial_point_no in range(3):
                for k in keys:
                    noise = noise_generator.get_noise(iter+1, trial_point_no, k)*1e4
                    fp.write("%.2e\n"%noise)

def get_noise_for_batch_test():
    noise_filename = "noise/noise1.txt"

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)
    for n in [0, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10]:
        noise_generator_new = NoiseGenerator()
        for k, v in noise_generator.history_noise.items():
            noise_generator_new.history_noise[k]=v*n
        noise_generator_new.save_noise("noise/batch_noise_%s.txt"%(str(n)))


if __name__ == "__main__":
    get_noise_for_batch_test()