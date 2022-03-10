from typing import Tuple, Generic, Callable, TypeVar, Any, List, Union
import re

CursorIdx = int
A = TypeVar("A")
B = TypeVar("B")


class ParseFailure(Exception):
    pass


class Parser(Generic[A]):
    parse: Callable[[str, CursorIdx], Tuple[A, CursorIdx]]

    def __init__(self, p: Callable[[str, CursorIdx], Tuple[A, CursorIdx]]):
        self.parse = p

    def run_parser(self, s: str) -> Union[Tuple[A, None], Tuple[None, ParseFailure]]:
        try:
            result = self.parse(s, 0)[0]
            return result, None
        except ParseFailure as e:
            return None, e


# This is the Functor of Parser that simply transform the result
# of any parse on the given parser by the given callable.
def parser_map(a: "Parser[A]", f: Callable[[A], B]) -> "Parser[B]":
    def map_result(s: str, c: CursorIdx):
        v, c = a.parse(s, c)
        return f(v), c

    return Parser(map_result)


def parser_lift(v: A) -> "Parser[A]":
    return Parser(lambda s, c: (v, c))


def parser_apply(f: "Parser[Callable[[A], B]]", a: "Parser[A]") -> "Parser[B]":
    def apply_result(s: str, c: CursorIdx) -> Tuple[B, CursorIdx]:
        next_work, c = f.parse(s, c)
        next_a, c = a.parse(s, c)
        return next_work(next_a), c

    return Parser(apply_result)


def parser_bind(f: "Parser[Callable[[A], Parser[B]]]", a: "Parser[A]") -> "Parser[B]":
    def bind_result(s: str, c: CursorIdx) -> Tuple[B, CursorIdx]:
        next_work, c = f.parse(s, c)
        next_a, c = a.parse(s, c)
        next_parser: "Parser[B]" = next_work(next_a)
        return next_parser.parse(s, c)

    return Parser(bind_result)


def either_parser(a: Parser[A], *alternatives: Parser[A]) -> Parser[A]:
    def try_one_then(s: str, c: CursorIdx) -> Tuple[A, CursorIdx]:
        try:
            return a.parse(s, c)
        except ParseFailure as e:
            last_exception = e

        for alt in alternatives:
            try:
                return alt.parse(s, c)
            except ParseFailure as e:
                last_exception = e

        raise last_exception

    return Parser(try_one_then)


def parse_many(parser: Parser[A]) -> Parser[List[A]]:
    def do_parsing(input_text: str, cursor: CursorIdx):
        items: List[A] = []

        while True:
            try:
                result, new_cursor = parser.parse(input_text, cursor)
                if new_cursor <= cursor:
                    raise ParseFailure(
                        "Potential infinite recursion in parse_many, cursor did not move forward."
                    )
                cursor = new_cursor
                items.append(result)
            except ParseFailure as failure:
                break

        return items, cursor

    return Parser(do_parsing)


# "Lift" a parsing value into a cell (singleton array) if it succeeds,
# otherwise provide an empty array.
def parser_maybe(parser: Parser[A]) -> Parser[List[A]]:
    return either_parser(parser_map(parser, lambda v: [v]), parser_lift([]))


def parse_many_joined_by(parser: Parser[A], sep: Parser[Any]) -> Parser[List[A]]:
    def with_head(head: List[A]) -> Parser[List[A]]:
        return parser_map(parser_maybe(parser), lambda tail: head + tail)

    each_head: Parser[A] = take_left(parser, sep)
    head_parser: Parser[List[A]] = parse_many(each_head)
    return parser_bind(parser_lift(with_head), head_parser)


def regex_parser(pattern: re.Pattern) -> Parser[str]:
    def do_parsing(input_text: str, cursor: CursorIdx):
        m = pattern.match(input_text, cursor)
        if m is None:
            raise ParseFailure(f"No match found for {pattern} at position {cursor}")

        start, end = m.span()
        return input_text[start:end], end

    return Parser(do_parsing)


def lift_parse_failure(message: str) -> Parser[Any]:
    def do_parsing(input_text: str, cursor: CursorIdx):
        raise ParseFailure(message)

    return Parser(do_parsing)


def literal_parser(literal: str) -> Parser[str]:
    def do_parsing(input_text: str, cursor: CursorIdx):
        next = input_text[cursor : cursor + len(literal)]
        if next != literal:
            raise ParseFailure(f"No match found for '{literal}' at position {cursor}")

        return literal, cursor + len(literal)

    return Parser(do_parsing)


def take_right(l: Parser[Any], r: Parser[A]) -> Parser[A]:
    def compose(_):
        def get_value(v):
            return v

        return get_value

    return parser_apply(parser_apply(parser_lift(compose), l), r)


def take_left(l: Parser[A], r: Parser[Any]) -> Parser[A]:
    def compose(v):
        def skip(_):
            return v

        return skip

    return parser_apply(parser_apply(parser_lift(compose), l), r)


def surrounded_by(l: Parser[Any], p: Parser[A], r: Parser[Any]) -> Parser[A]:
    return take_left(take_right(l, p), r)


whitespace_skipper: Parser[str] = regex_parser(re.compile(r"\s*"))


def token(p: Parser[A], bounds: Parser[Any] = whitespace_skipper) -> Parser[A]:
    return surrounded_by(bounds, p, bounds)


int_parser: Parser[int] = parser_map(regex_parser(re.compile(r"[-+]?[0-9]+")), int)
float_parser: Parser[float] = parser_map(regex_parser(re.compile(r"0\.\d\d+")), float)
eof_parser: Parser[str] = regex_parser(re.compile(r"$"))
