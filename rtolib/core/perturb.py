import copy
from rtolib.core import ProblemDescription

__name__=['SimpleFiniteDiffPerturbation']

class PerturbationMethod():
    def __init__(self):
        self.number_of_data_points = 0

    def get_trial_points(self, base_point):
        return NotImplementedError()

    def calculate_modifiers(self, plant_output_data, model_output_data):
        return NotImplementedError()


class SimpleFiniteDiffPerturbation(PerturbationMethod):
    def __init__(self, ffd_step, problem_description):
        '''
        :param self.ffd_step: dict.
         Stepsizes can be negative. the sign indicates the preferred direction
        :return:
        '''
        assert isinstance(problem_description, ProblemDescription)
        self.ffd_step = ffd_step
        self.mv_bounds = {}
        self.mvs = problem_description.symbol_list['MV']
        for v in problem_description.symbol_list['MV']:
            self.mv_bounds[v]=problem_description.bounds[v]
        self.number_of_data_points = len(self.mvs)+1

    def get_trial_points(self, base_point):
        '''

        :param base_point: dict
        :return: list of dicts
        '''
        ret = []
        ret.append(base_point)
        for k in self.mvs:
            point = copy.deepcopy(base_point)
            if self.ffd_step[k] > 0:
                if self.mv_bounds[k][1] is not None and point[k]+self.ffd_step[k] <= self.mv_bounds[k][1]:
                    point[k]+=self.ffd_step[k]
                elif self.mv_bounds[k][0] is not None and point[k] - self.ffd_step[k] >= self.mv_bounds[k][0]:
                    point[k] -= self.ffd_step[k]
                elif self.mv_bounds[k][1] is None and self.mv_bounds[k][0] is None:
                    point[k] += self.ffd_step[k]
                else:
                    raise Exception("The self.mv_bounds for %s is too tight or the self.ffd_step is too big."%k)
            elif self.ffd_step[k] <0:
                if self.mv_bounds[k][0] is not None and point[k]+self.ffd_step[k] >= self.mv_bounds[k][0]:
                    point[k]+=self.ffd_step[k]
                elif self.mv_bounds[k][1] is not None and point[k] - self.ffd_step[k] <= self.mv_bounds[k][1]:
                    point[k] -= self.ffd_step[k]
                elif self.mv_bounds[k][1] is None and self.mv_bounds[k][0] is None:
                    point[k] += self.ffd_step[k]
                else:
                    raise Exception("The self.mv_bounds for %s is too tight or the self.ffd_step is too big."%k)
            else:
                raise ValueError("FFD step has zero element.")
            ret.append(point)
        return ret

    def calculate_modifiers(self, plant_output_data, model_output_data):
        '''

        :param plant_output_data: list of dicts
        :param model_output_data: list of dicts
        :return: dict
        '''
        ret = {}
        cvs = model_output_data[0].keys()
        for cv in cvs:
            ret[(cv,None)] = plant_output_data[0][cv] - model_output_data[0][cv]
            for idx,mv in enumerate(self.mvs):
                ret[(cv, mv)] = (plant_output_data[idx+1][cv] - plant_output_data[0][cv])/self.ffd_step[mv]- \
                                (model_output_data[idx + 1][cv] - model_output_data[0][cv]) / self.ffd_step[mv]
        return ret