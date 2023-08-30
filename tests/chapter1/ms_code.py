
class InitializeSimulation:
    def __init__(self,
                 *args,
                 **kwargs,
                 ):
        super().__init__(*args, **kwargs)
        
        print("Initialize Simulation")


class Outputs:
    def __init__(self,
                 *args,
                 **kwargs):
        super().__init__(*args, **kwargs)

        print("Outputs")


class Utilities:
    def __init__(self,
                *args,
                **kwargs):
        super().__init__(*args, **kwargs)

        print("Utilities")


class MolecularDynamics(InitializeSimulation, Utilities, Outputs):
    def __init__(self,
                *args,
                **kwargs,
                ):
        super().__init__(*args, **kwargs)

        print("Start molecular dynamics simulation")


class MonteCarlo(InitializeSimulation, Utilities, Outputs):
    def __init__(self,
                 *args,
                 **kwargs,
                 ):
        super().__init__(*args, **kwargs)

        print("Start Monte Carlo simulation")