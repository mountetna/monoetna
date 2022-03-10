from typing import Any, Dict, Iterable, Tuple, Union
from serde import to_dict


def encode_as_multipart(v: Any) -> Dict:
    """
    Can be passed as the files= parameter of a requests.Session post call, for instance.
    """
    result = {}

    if not isinstance(v, dict):
        if hasattr(v, 'to_dict'):
            v = v.to_dict()
        else:
            v = to_dict(v)

    for k, v in _encode_as_multipart(v, "", True):
        result[k] = v

    return result


def _encode_as_multipart(
    v: Any, base_key: str, is_root: bool = False
) -> Iterable[Tuple[str, Union[str, Tuple[str, bytes, str]]]]:
    if isinstance(v, dict):
        for k, v in v.items():
            yield from _encode_as_multipart(v, k if is_root else f"{base_key}[{k}]")
    elif isinstance(v, list):
        for i, item in enumerate(v):
            if isinstance(v, dict):
                # This is necessary to ensure that arrays of hashes that have hetergenous keys still get parsed correctly
                # Since the only way to indicate a new entry in the array of hashes is by re-using a key that existed in
                # the previous hash.
                yield from _encode_as_multipart(i, f"{base_key}[][_idx]")

            yield from _encode_as_multipart(v, f"{base_key}[]")
    else:
        if not base_key:
            raise ValueError("base_key cannot be empty for scalars!")

        if isinstance(v, bytes):
            yield base_key, ("blob", v, "application/octet-stream")
        else:
            yield base_key, (None, str(v))
