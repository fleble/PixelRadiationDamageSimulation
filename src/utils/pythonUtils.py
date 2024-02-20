from importlib import import_module


def remove_double_characters(string, character):
    while 2*character in string:
        string = string.replace(2*character, character)
    return string


def import_module_from(file_name):
    """To do "import moduleName as module" in python script.

    Use as:
    module = import_module_from(python_file_name)

    Args:
        file_name (str)
    """

    file_name = remove_double_characters(file_name, "/")
    module_name = file_name.replace(".py", "").replace("/", ".")
    module = import_module(module_name)

    return module

