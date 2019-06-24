"""Defines a TaskSpec object"""
from pathlib import Path


class TaskSpec:
    """
    @brief      Class used to define a task spec in a standardized way
    """

    def __init__(
        self,
        func,
        func_fw_out,
        task_with_atoms_obj,
        func_kwargs=None,
        func_fw_out_kwargs=None,
        args=None,
        inputs=None,
        make_abs_path=False,
    ):
        """
        TaskSpec Constructor
        Args:
            func (str or Function): Function to be wrapped into a PyTask
            func_fw_out (str or Function): Function that converts the outputs of func into FWActions
            task_with_atoms_obj (bool): True if calculating using an ASE Atoms object
            func_kwargs (dict): kwargs for func
            func_fw_out_kwargs: kwargs for func_fw_out
            args (list): args for func
            inputs (list): args for func stored in the FireWorks database
                           (spec keys that are append to args)
            make_abs_path (bool): If True make all paths absolute
        """
        if not isinstance(func, str):
            func = f"{func.__module__}.{func.__name__}"
        if not isinstance(func_fw_out, str):
            func_fw_out = f"{func_fw_out.__module__}.{func_fw_out.__name__}"
        self.func = func
        self.func_fw_out = func_fw_out
        self.task_with_atoms_obj = task_with_atoms_obj

        if func_kwargs:
            if "workdir" in func_kwargs and make_abs_path:
                func_kwargs["workdir"] = str(Path(func_kwargs["workdir"]).absolute())
            self.func_kwargs = func_kwargs
        else:
            self.func_kwargs = {}

        if func_fw_out_kwargs:
            if "workdir" in func_fw_out_kwargs and make_abs_path:
                func_fw_out_kwargs["workdir"] = str(
                    Path(func_fw_out_kwargs["workdir"]).absolute()
                )
            self.func_fw_out_kwargs = func_fw_out_kwargs
        else:
            self.func_fw_out_kwargs = dict()

        if args:
            self.args = args
        else:
            self.args = list()

        if inputs:
            self.inputs = inputs
        else:
            self.inputs = list()

    def get_pt_args(self):
        """ get the PyTask args for the task """
        if self.task_with_atoms_obj:
            return [
                self.func,
                self.func_fw_out,
                self.func_kwargs,
                self.func_fw_out_kwargs,
                *self.args,
            ]
        return [self.func, self.func_fw_out, *self.args]

    def get_pt_kwargs(self, fw_settings):
        """ get the PyTask kwargs for the task """
        if not fw_settings:
            fw_settings = {}

        if self.task_with_atoms_obj:
            return {"fw_settings": fw_settings}

        to_ret = dict(self.func_kwargs, **self.func_fw_out_kwargs)
        to_ret["fw_settings"] = fw_settings

        return to_ret

    def get_pt_inputs(self):
        """ get the PyTask inputs for the task """
        return self.inputs
