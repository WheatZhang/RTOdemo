

class BlackBoxModel(object):
    def build(self):
        return NotImplementedError

    def simulate(self, input_dict):
        return NotImplementedError


class BlackBoxModelWithModifiers(BlackBoxModel):
    def __init__(self, black_box_model, mvs, cvs):
        '''

        :param black_box_model: BlackBoxModel
        :param mvs: list of strings
        :param cvs: list of strings
        :return:
        '''
        assert isinstance(black_box_model, BlackBoxModel)
        self.base_bb_model = black_box_model
        self.modifier_mvs = mvs
        self.modifier_cvs = cvs
        self.modifiers_value = {}
        for name in self.modifier_names_iterator():
            self.modifiers_value[name] = 0
        self.base_input = {}

    def build(self):
        self.base_bb_model.build()

    def modifier_names_iterator(self):
        for cv in self.modifier_cvs:
            yield cv + "_eps"
            for mv in self.modifier_mvs:
                yield mv+"_"+cv + "_lam"

    def set_base_point(self, base_point):
        self.base_input = base_point

    def set_modifiers(self, modifiers_value):
        for name in modifiers_value.keys():
            if name in self.modifiers_value.keys():
                self.modifiers_value[name] = modifiers_value[name]
            else:
                raise KeyError("unknown modifier name %s"%name)
    def simulate(self, input_dict):
        output_dict = self.base_bb_model.simulate(input_dict)
        for cv in self.modifier_cvs:
            output_dict += self.modifiers_value[cv + "_eps"]
            temp = 0
            for mv in self.modifier_mvs:
                temp += self.modifiers_value[mv+"_"+cv + "_lam"] * (input_dict[mv]-self.base_input[mv])
            output_dict += temp
        return output_dict
