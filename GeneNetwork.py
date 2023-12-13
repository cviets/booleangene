#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 13:15:58 2023

@author: chrisviets
"""

import doctest
import random
import csv


def remove_extra_parens(lst):
    """
    Given a list representing a formula, removes unneeded parentheses around
    literals (eg, ( NOT A ) AND ( B ) can be represented as NOT A AND B )

    Input: list
    Output: None

    >>> lst = ['(', 'A', ')']
    >>> remove_extra_parens(lst)
    >>> print(lst)
    ['A']

    >>> lst = ['(', 'NOT', 'A', ')']
    >>> remove_extra_parens(lst)
    >>> print(lst)
    ['NOT', 'A']

    """
    remove = True
    while remove:
        remove = set()
        for i, elt in enumerate(lst):
            if elt == "(":
                end_idx = find_closing_parens(lst, i)
                smaller_lst = lst[i + 1 : end_idx]
                if not ("AND" in smaller_lst or "OR" in smaller_lst):
                    remove.update({i, end_idx})
        for idx in sorted(remove, reverse=True):
            del lst[idx]


def find_closing_parens(lst, idx):
    """
    Given a list representing a formula and an index in the list (where the element
    at 'idx' must be an open parenthesis), returns the index in the list corresponding
    to the closing parenthesis. If none found, returns -1

    Input: list, int (index in list)
    Output: int (index in list)
    """
    assert lst[idx] == "("
    count = 0
    for i in range(idx, len(lst)):
        elt = lst[i]
        if elt == "(":
            count += 1
        if elt == ")":
            count -= 1
        if count == 0:
            return i
    return -1


def parse(formula, parameters):
    """
    Evaluates formula with variable values given by parameters dictionary.
    Note that formula must not contain special Python characters (e.g. '/', ',',
    '.', '+'); only (variable names composed of letters, numbers, and '_'s) (', ')',
    'AND', 'OR', 'NOT' and spaces as separators are allowed.

    Input: formula (list), parameters (dict)
    Output: bool (truth value of formula given parameters)
    """

    newf = formula.replace("AND", "and")
    newf = newf.replace("OR", "or")
    newf = newf.replace("NOT", "not")

    out = eval(newf, parameters | {"__builtins__": __builtins__})
    return out


def read_expressions_from_txt(path_to_txt_file):
    """
    Returns dictionary that maps genes to their corresponding expressions in input
    text file.
    Makes the following replacements: {',': '_', '/': '_', '+', '__', '-': ''}

    Input: string (path to txt file)
    Output: dict (keys = gene names, values = string expressions)

    """
    out = {}
    reader = open(path_to_txt_file, "r")
    for line in reader:
        if line[-1] == "\n":
            line = line[:-1]
        var, formula = line.split(" = ")
        out[
            var.replace(",", "_").replace("/", "_").replace("+", "__").replace("-", "")
        ] = (
            formula.replace(",", "_")
            .replace("-", "")
            .replace("+", "__")
            .replace("/", "_")
        )
    return out


def read_external_from_txt(path_to_txt_file):
    """
    Returns set containing external signals given text file.
    Makes the following replacements: {',': '_', '/': '_', '+', '__', '-': ''}

    Input: string (path to file)
    Output: set
    """
    out = set()
    reader = open(path_to_txt_file, "r")
    for line in reader:
        if line[-1] == "\n":
            line = line[:-1]
        out.add(
            line.replace(",", "_").replace("/", "_").replace("-", "").replace("+", "__")
        )
    return out


def split_formula(formula):
    """
    Given a string formula, represents formula as a 2-element tuple. 0th element
    is a 2-element tuple containing formula on left and right side of operator.
    1st element is 'AND' or 'OR', representing the main operator of the formula.
    Applies DeMorgan's rules on formula first so main operator is always 'OR' or
    'AND'.

    'link' is the inverse function of 'split_formula'

    Input: string
    Output: tuple (may be nested)

    >>> split_formula("A OR B")
    (('A', 'B'), 'OR')

    >>> split_formula("NOT A AND NOT B")
    (('NOT A', 'NOT B'), 'AND')

    >>> split_formula("( A OR B ) AND NOT C")
    (((('A', 'B'), 'OR'), 'NOT C'), 'AND')

    >>> split_formula("A")
    'A'

    >>> split_formula("NOT ( Am OR A )")
    (('NOT Am', 'NOT A'), 'AND')

    >>> split_formula("( Am OR A OR R )")
    (('Am', (('A', 'R'), 'OR')), 'OR')

    >>> split_formula("NOT ( A OR R )")
    (('NOT A', 'NOT R'), 'AND')

    >>> split_formula("NOT ( Am OR A OR R )")
    (('NOT Am', (('NOT A', 'NOT R'), 'AND')), 'AND')

    >>> split_formula("NOT ( NOT Am OR A )")
    (('Am', 'NOT A'), 'AND')

    """

    formula_lst = [elt for elt in formula.split(" ") if elt != ""]
    remove_extra_parens(formula_lst)
    remove_double_nots(formula_lst)

    if (
        formula_lst[0] == "("
        and find_closing_parens(formula_lst, 0) == len(formula_lst) - 1
    ):
        return split_formula(" ".join(formula_lst[1:-1]))

    if len(formula_lst) < 3:
        return " ".join(formula_lst)

    if (
        formula_lst[0] == "NOT"
        and formula_lst[1] == "("
        and find_closing_parens(formula_lst, 1) == len(formula_lst) - 1
    ):
        count = 0
        modified_lst = formula_lst[2:-1]

        for i, elt in enumerate(modified_lst):
            if elt == "(":
                count += 1
            if elt == ")":
                count -= 1
            if count == 0 and (elt == "AND" or elt == "OR"):
                subformula_1 = "NOT ( " + " ".join(modified_lst[:i]) + " )"
                not_operator = elt
                subformula_2 = "NOT ( " + " ".join(modified_lst[i + 1 :]) + " )"

                if not_operator == "AND":
                    operator = "OR"
                else:
                    operator = "AND"

                return (
                    (split_formula(subformula_1), split_formula(subformula_2)),
                    operator,
                )

    count = 0

    for i, elt in enumerate(formula_lst):
        if elt == "(":
            count += 1
            continue
        if elt == ")":
            count -= 1
            continue
        if count == 0 and (elt == "AND" or elt == "OR"):
            subformula_1 = formula_lst[:i]
            operator = elt
            subformula_2 = formula_lst[i + 1 :]

            return (
                (
                    split_formula(" ".join(subformula_1)),
                    split_formula(" ".join(subformula_2)),
                ),
                operator,
            )


def reduce_formula(branched_formula, variable):
    """
    Given formula as tuple and variable, returns new tuple formula assuming 
    variable is True. Used as helper for 'reduce' function.
    
    >>> reduce_formula((("True", "NOT C"), "OR"), "A")
    'True'

    >>> reduce_formula((("False", "NOT C"), "AND"), "A")
    'False'

    >>> reduce_formula((("A", "B"), "OR"), "A")
    'True'

    >>> reduce_formula((("NOT A", "B"), "AND"), "A")
    'False'

    >>> reduce_formula("C", "A")
    'C'

    >>> reduce_formula((((('A', 'B'), 'OR'), 'NOT C'), 'AND'), "A")
    'NOT C'

    >>> reduce_formula((('A', 'R'), 'OR'), 'C')
    (('A', 'R'), 'OR')

    >>> reduce_formula((('Am', (('A', 'R'), 'OR')), 'OR'), "C")
    (('Am', (('A', 'R'), 'OR')), 'OR')

    >>> reduce_formula((((('C', 'NOT R'), 'AND'), 'NOT Rm'), 'AND'), 'C')
    (('NOT R', 'NOT Rm'), 'AND')

    """

    if isinstance(branched_formula, str):
        if branched_formula == variable:
            return "True"
        if branched_formula == "NOT " + variable:
            return "False"
        return branched_formula

    operator = branched_formula[1]
    subformula_1 = branched_formula[0][0]
    subformula_2 = branched_formula[0][1]

    if "True" in branched_formula[0] or variable in branched_formula[0]:
        if operator == "OR":
            return "True"
        if operator == "AND":
            if subformula_1 in {"True", variable}:
                return reduce_formula(subformula_2, variable)
            if subformula_2 in {"True", variable}:
                return reduce_formula(subformula_1, variable)

    if "False" in branched_formula[0] or "NOT " + variable in branched_formula[0]:
        if operator == "AND":
            return "False"
        if operator == "OR":
            if subformula_1 in {"False", "NOT " + variable}:
                return reduce_formula(subformula_2, variable)
            if subformula_2 in {"False", "NOT " + variable}:
                return reduce_formula(subformula_1, variable)

    if isinstance(subformula_1, str) and isinstance(subformula_2, str):
        return branched_formula

    out = (
        (
            reduce_formula(subformula_1, variable),
            reduce_formula(subformula_2, variable),
        ),
        operator,
    )
    if out == branched_formula:
        return out
    return reduce_formula(out, variable)


def link(branched_formula, outermost=True):
    """
    Inverse of split_formula. Given a tuple representing a formula, returns formula
    in string format.

    Input: tuple (may be nested) (outermost is optional argument only needed
                                  for recursion)
    Output: string

    >>> link("A")
    'A'

    >>> link((('A', 'B'), 'OR'))
    '( A OR B )'

    >>> link((('NOT Am', (('NOT A', 'NOT R'), 'AND')), 'AND'))
    'NOT Am AND ( NOT A AND NOT R )'

    >>> link((('C', (('NOT A', 'B'), 'AND')), 'OR'))
    'C OR ( NOT A AND B )'

    >>> link((('R', (('NOT Am', (('NOT A', 'NOT R'), 'AND')), 'AND')), 'OR'))
    'R OR ( NOT Am AND ( NOT A AND NOT R ) )'

    """
    if isinstance(branched_formula, str):
        return branched_formula

    operator = branched_formula[1]
    subformula_1 = branched_formula[0][0]
    subformula_2 = branched_formula[0][1]

    if isinstance(subformula_1, str) and isinstance(subformula_2, str):
        return "( " + subformula_1 + " " + operator + " " + subformula_2 + " )"
    if not outermost:
        return (
            "( "
            + link(subformula_1, outermost=False)
            + " "
            + operator
            + " "
            + link(subformula_2, outermost=False)
            + " )"
        )
    else:
        return (
            link(subformula_1, outermost=False)
            + " "
            + operator
            + " "
            + link(subformula_2, outermost=False)
        )


def falsify_variable(formula, variable):
    """
    Returns modified formula such that 'variable' in formula is made False at
    every instance.
    
    Input: formula (string), variable (string)
    Output: new formula (string)
    
    >>> falsify_variable("NOT D", "D")
    'D'

    >>> falsify_variable("( D AND NOT C ) OR ( A OR B OR C )", "C")
    '( D AND C ) OR ( A OR B OR NOT C )'

    >>> falsify_variable("A AND ( B OR C OR NOT A )", "A")
    'NOT A AND ( B OR C OR A )'
    """

    formula_lst = [elt for elt in formula.split(" ") if elt != ""]
    remove_extra_parens(formula_lst)

    out = [elt if elt != variable else "NOT " + elt for elt in formula_lst]
    out = " ".join(out)
    out = [elt for elt in out.split(" ") if elt != ""]
    for i in range(len(out) - 2, -1, -1):
        elt = out[i]
        if elt == "NOT" and out[i + 1] == "NOT":
            del out[i + 1]
            del out[i]

    return " ".join(out)


def remove_double_nots(lst):
    """
    Given list representing formula, removes pairs of consecutive 'NOT's.
    Mutates original list, returns None.
    """

    for i in range(len(lst) - 2, -1, -1):
        elt = lst[i]
        if elt == "NOT" and lst[i + 1] == "NOT":
            del lst[i + 1]
            del lst[i]


def reduce(formula, variable, falsify=False):
    """
    Given formula and variable, returns modified formula assuming variable is True
    or assuming variable is False (if falsify is True).
    
    >>> Rm = "( ( R ) AND NOT ( ( Am )  OR ( A ) ) ) OR  NOT ( Am OR A OR R ) "
    >>> reduce(Rm, "A")
    'False'

    >>> reduce("A AND B", 'A')
    'B'

    >>> reduce("( A OR B OR C )", 'D')
    'A OR ( B OR C )'

    >>> reduce("( ( ( C  ) AND NOT ( R  )  ) AND NOT ( Rm  ) ) ", "A")
    '( C AND NOT R ) AND NOT Rm'

    >>> reduce("( ( ( C  ) AND NOT ( R  )  ) AND NOT ( Rm  ) ) ", "C")
    'NOT R AND NOT Rm'

    ADD MORE DOCTESTS HERE. THIS FUNCTION IS IMPORTANT
    """

    if falsify:
        return reduce(falsify_variable(formula, variable), variable, False)

    branched = split_formula(formula)
    out = reduce_formula(branched, variable)

    if isinstance(out, str):
        return out

    out = link(out)
    formula_lst = [elt for elt in out.split(" ") if elt != ""]
    remove_extra_parens(formula_lst)
    if (
        formula_lst[0] == "("
        and find_closing_parens(formula_lst, 0) == len(formula_lst) - 1
    ):
        formula_lst = formula_lst[1:-1]
    return " ".join(formula_lst)


def get_inputs(formula):
    return tuple(
        elt
        for elt in formula.split(" ")
        if elt not in {"", "NOT", "OR", "AND", "(", ")"}
    )


def get_input_assignments(literals):
    input_as_string = "0b" + "0" * len(literals)
    while input_as_string != "0b1" + "0" * len(literals):
        yield {
            literals[i]: bool(int(char))
            for i, char in enumerate(input_as_string[2:])
            if literals[i] != "__builtins__"
        }
        input_as_string = bin(int(input_as_string, 2) + 1)
        while len(input_as_string) < len(literals) + 2:
            input_as_string = "0b0" + input_as_string[2:]


def translate_state_to_string(state):
    out = ""
    for elt in state.values():
        out = out + str(int(elt))
    return out


def get_truth_table(formula="", read_from_file=False):
    """
    Given either string formula xor csv file containing truth table, returns
    truth table as dict that maps bitstrings to truth, and returns tuple of gene
    names (strings) to indicate order of bitstrings.
    
    Input: formula (string) XOR csv file path (string)
    Output: truth table (dict) and variable order (tuple)
    """
    out = {}

    if read_from_file:
        with open(read_from_file, newline="") as csvfile:
            heading = (
                next(csvfile).replace("/", "_").replace("-", "").replace("+", "__")
            )
            variables = (
                " ".join(heading.split(" ")[:-1]).replace(", ", " ").replace(",", "_")
            )
            if variables[-1] == "_":
                variables = variables[:-1]
            variables = tuple(variables.split(" "))
            reader = csv.reader(csvfile, delimiter=",")
            for line in reader:
                in_out = "".join(line)
                in_out = in_out.replace(" ", "")
                inp = in_out[:-1]
                outp = False if in_out[-1] == "0" else True
                out[inp] = outp

    else:
        literals = get_inputs(formula)
        for tva in get_input_assignments(literals):
            state_string = translate_state_to_string(tva)
            out[state_string] = parse(formula, tva)

        variables = tuple(elt for elt in tva)
    return out, variables


def check_canalizability(truth_table, idx, value):
    """
    Helper function for 'canalize'. Checks truth table (dict) to see if gene at 
    'idx' in keys with truth-value given by 'value' is canalizing. If so, returns
    truth value that formula takes on when canalized. If not, returns None
    
    Input: truth_table (dict), idx (int), value (bool)
    Output: bool if canalizible, else None
    """
    value = int(value)
    smaller_table = {
        elt: val for elt, val in truth_table.items() if int(elt[idx]) == value
    }

    permitted_vals = set(smaller_table.values())
    if len(permitted_vals) == 1:
        return tuple(permitted_vals)[0]
    return None


def assign_single_variable(num_total):
    """
    Helper function for canalize. Gemnerator that yields gene number (0 to num_total - 1) 
    and truth value (0 or 1), which are then fed into check_canalizability

    """
    count = 0

    while count < num_total * 2:
        yield count // 2, count % 2
        count += 1


def canalize(formula, read_tt_from_csv=False):
    """
    Returns formula's canalizing variable and its canalizing truth value (eg, "A"
    or "NOT A"), as well as value that 'formula' takes when the canalizing condition
    is met. For instance, ("NOT A", True) indicates that 'formula' is True when
    A is False. Function finds canalizing values by first computing formula's 
    truth table. This step can be skipped if user inputs truth table csv file.
    
    Input: formula (string), read_tt_from_csv=False (optional, allows user to 
                                                     input truth table as csv)
    Output: canalizing variable (string) and formula value when canalized (bool)
    """
    def extra(num):
        if num == 0:
            return "NOT "
        return ""

    truth_table, variables = get_truth_table(
        formula=formula, read_from_file=read_tt_from_csv
    )
    for idx, assignment in assign_single_variable(len(variables)):
        val = check_canalizability(truth_table, idx, assignment)
        if val is not None:
            return extra(assignment) + variables[idx], val


def nested_canalize(formula, read_tt_from_csv=False):
    """
    NOTE: truth table is recomputed for each level in hierarchy. To save time,
    future versions of this function should reuse truth table from top of hierarchy
    instead of recomputing each time.
    
    Returns tuple describing hierarchy of canalization where each element of tuple
    corresponds to 2-tuple of form output by 'canalize' function. Also returns
    subformula when all conditions in hierarchy are not met.
    
    For instance, the output (('A', False), ('NOT B', True)), 'True' indicates
    that if 'A' is True, then formula is False. Then, assuming 'A' is False, 
    if 'NOT B' is True then formula is True. If 'A' and 'NOT B' are both False,
    then the formula is True
    
    Input: formula (string), path to csv file (optional, string)
    Output: tuple of tuples, string (subformula)

    """
    new_formula = formula[:]
    out = tuple()
    canal = canalize(new_formula, read_tt_from_csv=read_tt_from_csv)

    while canal is not None and new_formula not in {"True", "False"}:
        out += (canal,)
        if canal[0][:4] == "NOT ":
            new_formula = reduce(new_formula, canal[0][4:])
        else:
            new_formula = reduce(new_formula, canal[0], falsify=True)
        canal = canalize(new_formula)

    return out, new_formula


class BooleanNetwork:
    def __init__(self, expressions, external):
        # self.state and self.expressions are dicts
        # (map var to truth value or var to formula)

        self.state = {gene: None for gene in expressions}
        self.expressions = expressions
        self.external = external

    def initialize_external(self, ext_dict):
        self.external_values = ext_dict

    def initial_condition(self, state):
        self.state = state

    def step(self):
        new_state = {var: truth for var, truth in self.state.items()}
        for key, val in self.expressions.items():
            new_state[key] = parse(val, self.state | self.external_values)
        self.state = new_state

    def find_attractors(self, stable_only=False, n=1000, t=100):
        """
        Returns network's attractors (including stable fixed points and oscillatory
        cycles). Excludes oscillatory cycles if stable_only is True. Finds attractors
        by running n simulations that each runs for t time steps.
        
        Input (all optional): stable_only (bool), n (int), t (int)
        Output: dict where keys are tuples corresponding to attractors and values
        are ints showing how many times in the n simulations that attractor was
        reached.

        """
        out = {}
        for __ in range(n):
            params = {}
            for gene in self.expressions:
                params[gene] = bool(random.randint(0, 1))

            self.initial_condition(params)
            history = [translate_state_to_string(self.state)]

            for _ in range(t):
                self.step()
                translated_state = translate_state_to_string(self.state)

                if translated_state in history:
                    idx = history.index(translated_state)
                    cycle = tuple(history[idx:])

                    if len(cycle) > 1 and stable_only:
                        break

                    found = False
                    for elt in out:
                        if set(elt) == set(cycle):
                            out[elt] += 1
                            found = True
                            break
                    if found:
                        break
                    out[cycle] = 1
                    break

                history.append(translated_state)
        return out


# class TruthTree:
#     def __init__(self, set_of_formulas):
#         self.formulae = set_of_formulas

#         if len(self.formulae) == 1:
#             formula = list(self.formulae)[0]
#             branches = split_formula(formula)

#             if isinstance(branches, str):
#                 self.children = None
#                 self.operator = None
#             else:
#                 self.operator = branches[1]
#                 self.subformula = branches[0]

#                 if self.operator == "OR":
#                     self.children = (
#                         TruthTree(self.subformula[0]),
#                         TruthTree(self.subformula[1]),
#                     )
#                 if self.operator == "AND":
#                     other = TruthTree({self.subformula[0], self.subformula[1]})
#                     self.children = (other,)

#         else:
#             for formula in self.formulae:
#                 branches = split_formula(formula)


if __name__ == "__main__":

    # doctest.testmod()

    # _doctest_flags = doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS
    # run all doctests
    # doctest.testmod(optionflags=_doctest_flags)

    # doctest.run_docstring_examples(
    #     link, globals(), optionflags=_doctest_flags, verbose=False
    # )

    external = read_external_from_txt(
        "/Users/chrisviets/Downloads/expr/external_components.ALL.txt"
    )
    expressions = read_expressions_from_txt(
        "/Users/chrisviets/Downloads/expr/expressions.ALL.txt"
    )

    net = BooleanNetwork(expressions, external)
    ext_dict = {}
    for elt in {'Le', 'Ge', 'Lem'}:
        ext_dict[elt] = bool(random.randint(0, 1))

    # # oscillations possible iff {'Ge': False, 'Lem': True, 'Le': False}
    ext_dict = {'Ge': False, 'Lem': True, 'Le': False}
    net.initialize_external(ext_dict)
    print(net.find_attractors(stable_only=False, n=1000, t=10))


    formula = '( NFKB1_MAP3K8_complex AND ( ( ( AKT2_phosphorylated ) ) )    )  OR ( IKK_complex AND ( ( ( NFKB1_MAP3K8_complex ) ) )    ) '
    print()
    print(nested_canalize(formula))
