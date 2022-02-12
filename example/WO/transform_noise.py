from rtolib.core import NoiseGenerator

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