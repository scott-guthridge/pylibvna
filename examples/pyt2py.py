#!/usr/bin/python3
from enum import Enum
import os
import re
import sys


class Parser:
    class State(Enum):
        DEFAULT = 0
        BRACKET = 1
        CURLY = 2

    def __init__(self):
        self.state = Parser.State.DEFAULT
        self.text = ""
        self.fstring = False
        self.indent = 0
        self.header = True
        self.outfile = None

    def _add(self, s):
        if self.state == Parser.State.BRACKET:
            print(s, end='')
        else:
            self.text += s

    def _flush(self):
        if self.text:
            if self.header:
                print('_file = None')
                self.header = False
            print('_code = ', end='')
            if self.fstring:
                print('f', end='')
            print('"""\\')
            print(self.text, end='')
            print('"""')
            print("print(_code, end='', file=_file)")
            self.text = ''
        self.fstring = False

    def parse(self, filename, file, tab=8):
        '''
        Parse a file
        '''
        line_number = 0
        for line in file:
            line = line.rstrip()
            cur = 0
            end = len(line)
            line_number += 1
            no_EOL = False

            # Expand tabs
            while (i := line.find('\t', cur)) != -1:
                cur = i
                npad = tab - (cur % tab)
                line = line[:cur] + (' ' * npad) + line[cur+1:]
                cur += npad - 1
                end += npad - 1
            cur = 0

            # Check and remove indentation of %[ block.
            if self.state == Parser.State.BRACKET and self.indent:
                if line != "" and line[:self.indent] != ' ' * self.indent:
                    raise IndentationError(f'{filename} (line {line_number}) '
                                           f'fatal: %[ block must be indented '
                                           f'from the opening')
                cur = self.indent

            # Parse line
            while cur < end:
                c = line[cur]

                # Advance to the next character
                def getc(noEnd=False):
                    nonlocal cur
                    nonlocal c
                    cur += 1
                    if cur != end:
                        c = line[cur]
                        return
                    if noEnd:
                        raise SyntaxError(f'{filename} (line {line_number}) '
                                          f'fatal: incomplete expression at '
                                          f'end of line')
                    c = '\n'

                # Backslash quote " and \ if going into self.text
                if c == '"' or c == '\\':
                    if self.state != Parser.State.BRACKET:
                        self._add('\\')
                    self._add(line[cur])
                    getc()
                    continue

                # Double {'s and }'s if we're in an f-string
                if c == '{' or c == '}':
                    if self.fstring:
                        self._add(c)
                    self._add(c)
                    getc()
                    continue

                # Handle %-escapes
                if c == '%':
                    percent_indent = cur
                    getc()
                    # %% is literal %
                    if c == '%':
                        self._add('%')
                        getc()
                        continue

                    # %# is comment
                    if c == '#':
                        if (self.state != Parser.State.DEFAULT
                                or line[:cur-1].strip() == ""):
                            no_EOL = True
                        # getc()
                        cur = end
                        continue

                    # %[ starts an embedded code block
                    if c == '[':
                        if self.state == Parser.State.BRACKET:
                            raise SyntaxError(f'{filename} '
                                              f'(line {line_number}) '
                                              'fatal: %[ blocks cannot '
                                              'be nested')

                        if self.state == Parser.State.CURLY:
                            raise SyntaxError(f'{filename} '
                                              f'(line {line_number}) '
                                              'fatal: expected %}')
                        self.indent = percent_indent
                        if percent_indent:
                            self.text = self.text[:-percent_indent]
                        self._flush()
                        # suppress _indent at top of output file
                        if not self.header:
                            print(f'_indent = {percent_indent}')
                        # getc()
                        cur = end             # skip rest of line
                        no_EOL = True
                        self.state = Parser.State.BRACKET
                        continue

                    # %] ends an embedded code block
                    if c == ']':
                        if self.state != Parser.State.BRACKET:
                            raise SyntaxError(f'{filename} '
                                              f'(line {line_number}) '
                                              'fatal: %] without %[')
                        # getc()
                        cur = end           # skip rest of line
                        no_EOL = True
                        if self.header:
                            print('_file = None')
                            self.header = False
                        self.state = Parser.State.DEFAULT
                        continue

                    # %{ starts an f-string {} expr
                    if c == '{':
                        if self.state == Parser.State.BRACKET:
                            raise SyntaxError(f'{filename} '
                                              f'(line {line_number}) '
                                              'fatal: %{ not allowed inside '
                                              '%[ ... %] block')

                        if self.state == Parser.State.CURLY:
                            raise SyntaxError(f'{filename} '
                                              f'(line {line_number}) '
                                              'fatal: %{ expressions '
                                              'cannot be nested')
                        if not self.fstring:
                            # Change all captured { and } to {{ and }}.
                            self.text = re.sub(r'([{}])', r'\1\1', self.text)
                            self.fstring = True
                        self._add('{')
                        getc()
                        self.state = Parser.State.CURLY
                        continue

                    # %} ends an f-string {} expr
                    if c == '}':
                        if self.state != Parser.State.CURLY:
                            raise SyntaxError(f'{filename} '
                                              f'(line {line_number}) '
                                              'fatal: %} without %{')
                        self._add('}')
                        getc()
                        self.state = Parser.State.DEFAULT
                        continue

                    # %O opens a new output file
                    if c == 'O':
                        if self.state != Parser.State.DEFAULT:
                            raise SyntaxError(f'{filename} '
                                              f'(line {line_number}) '
                                              'fatal: %O inside code block')
                        self._flush()
                        while True:
                            getc(noEnd=True)
                            if c != ' ':
                                break
                        if self.outfile:
                            print(f"_file.close()")
                        self.outfile = line[cur:]
                        print(f"_file = open('{self.outfile}', 'w')")
                        cur = end             # skip rest of line
                        no_EOL = True
                        self.header = False
                        continue

                    raise SyntaxError(f'{filename} '
                                      f'(line {line_number}) '
                                      f'unknown %-escape %{c}')

                # Else, a regular character
                self._add(c)
                getc()

            # Add end of line if not suppressed.
            if not no_EOL:
                self._add(os.linesep)

        # Check for EOF in a %[ block
        if self.state == Parser.State.BRACKET:
            raise SyntaxError(f'{filename} '
                              f'(line {line_number}) '
                              '%[ without matching %]')

        # Check for EOF in a %{ expression
        if self.state == Parser.State.CURLY:
            raise SyntaxError(f'{filename} '
                              f'(line {line_number}) '
                              '%{ without matching %}')

        # Flush any collected text
        self._flush()

        # Close the output file (if open)
        if self.outfile:
            print(f"_file.close()")
            self.outfile = None


parser = Parser()
try:
    if len(sys.argv) > 1:
        for filename in sys.argv[1:]:
            if filename == '-':
                parser.parse('<stdin>', sys.stdin)
            else:
                with open(filename, 'r', encoding='utf-8') as file:
                    parser.parse(filename, file)
    else:
        parser.parse('<stdin>', sys.stdin)
except Exception as e:
    print(e)
