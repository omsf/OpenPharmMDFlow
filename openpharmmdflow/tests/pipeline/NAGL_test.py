from openff.utilities.testing import skip_if_missing


@skip_if_missing("openff.toolkit.utils.nagl_wrapper")
def nagl_test():
    from openff.nagl_models import get_models_by_type
    from openff.nagl_models import validate_nagl_model_path
    from openff.toolkit import Molecule
    from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper

    nagl_model = get_models_by_type(model_type="am1bcc", production_only=True)[-1]
    model_path = validate_nagl_model_path(nagl_model)
    my_mol = Molecule.from_file("../polymer/HPMCAS-H-Trimer.sdf")
    my_mol.assign_partial_charges(
        toolkit_registry=NAGLToolkitWrapper(),
        partial_charge_method=model_path,
    )
