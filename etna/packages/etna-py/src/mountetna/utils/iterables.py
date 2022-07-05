from typing import Iterable, TypeVar, List

T = TypeVar("T")


def batch_iterable(iter: Iterable[T], batch_size: int) -> Iterable[Iterable[T]]:
    cur_batch: List[T] = []
    for item in iter:
        cur_batch.append(item)
        if len(cur_batch) >= batch_size:
            yield cur_batch
            cur_batch = []
    if cur_batch:
        yield cur_batch
