import dataclasses


# Modifies an existing dataclass instance by deconstructing any present values which are dictionaries but
# could be pact into another dataclass, enabling nested dataclasses.
import typing


def _unwrap_inner(instance):
    for f in dataclasses.fields(type(instance)):
        value = getattr(instance, f.name)

        if isinstance(value, dict):
            if not dataclasses.is_dataclass(f.type):
                continue
            new_value = f.type(**value)
            setattr(instance, f.name, new_value)
        elif isinstance(value, list):
            if typing.get_origin(f.type) == list:
                args = typing.get_args(f.type)
                if args and dataclasses.is_dataclass(args[0]):
                    inner_dc = args[0]
                    inner_fields = {f.name for f in dataclasses.fields(inner_dc)}
                    new_value = [(inner_dc(**({ k: v for k, v in row.items() if k in inner_fields })) if isinstance(row, dict) else row) for row in value]
                    setattr(instance, f.name, new_value)

class DataclassWithNesting:
    """
    Subclass this, in addition to using the @dataclasses.dataclass decorator to enable
    deserializing nested objects.
    """
    def __post_init__(self, *args, **kwds):
        _unwrap_inner(self)