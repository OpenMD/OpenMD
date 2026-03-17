"""
OMDFileVisitor.py

A concrete visitor for the OMD ANTLR4 grammar that builds a structured
Python representation of an OpenMD file.

Usage:
    from antlr4 import CommonTokenStream, FileStream
    from OMDLexer import OMDLexer
    from OMDParser import OMDParser
    from OMDVisitor import OMDVisitor

    input_stream = FileStream("my_file.omd")
    lexer = OMDLexer(input_stream)
    stream = CommonTokenStream(lexer)
    parser = OMDParser(stream)
    tree = parser.omdfile()

    visitor = OMDFileVisitor()
    result = visitor.visitOmdfile(tree)
"""

from OMDVisitor import OMDVisitor

# If ANTLR generated OMDParser, import context types from it.
# These imports assume antlr4-tools or antlr4 python runtime is installed
# and the grammar has been compiled with:
#   antlr4 -Dlanguage=Python3 -visitor OMD.g4
try:
    from OMDParser import OMDParser
except ImportError:
    OMDParser = None


class OMDFileVisitor(OMDVisitor):
    """
    Concrete visitor that walks the OMD parse tree and returns structured
    Python dicts/lists representing the file contents.

    The top-level visit(tree) call returns a dict of the form:
    {
        'assignments': { 'key': value, ... },
        'components':  [ { 'assignments': {...} }, ... ],
        'molecules':   [ { 'assignments': {...}, 'atoms': [...], ... }, ... ],
        'fragments':   [ ... ],
        'zconstraints': [ ... ],
        'restraints':  [ ... ],
        'flucq':       [ ... ],
        'rnemd':       [ ... ],
        'light':       [ ... ],
        'minimizer':   [ ... ],
    }
    """

    # -------------------------------------------------------------------------
    # Top level
    # -------------------------------------------------------------------------

    def visitOmdfile(self, ctx):
        result = {
            'assignments':  {},
            'components':   [],
            'molecules':    [],
            'fragments':    [],
            'zconstraints': [],
            'restraints':   [],
            'flucq':        [],
            'rnemd':        [],
            'light':        [],
            'minimizer':    [],
        }
        for stmt in ctx.statement():
            self._dispatch_statement(stmt, result)
        return result

    def _dispatch_statement(self, ctx, result):
        if ctx.assignment():
            key, value = self.visitAssignment(ctx.assignment())
            result['assignments'][key] = value
        elif ctx.componentblock():
            result['components'].append(self.visitComponentblock(ctx.componentblock()))
        elif ctx.moleculeblock():
            result['molecules'].append(self.visitMoleculeblock(ctx.moleculeblock()))
        elif ctx.fragmentblock():
            result['fragments'].append(self.visitFragmentblock(ctx.fragmentblock()))
        elif ctx.zconstraintblock():
            result['zconstraints'].append(self.visitZconstraintblock(ctx.zconstraintblock()))
        elif ctx.restraintblock():
            result['restraints'].append(self.visitRestraintblock(ctx.restraintblock()))
        elif ctx.flucqblock():
            result['flucq'].append(self.visitFlucqblock(ctx.flucqblock()))
        elif ctx.rnemdblock():
            result['rnemd'].append(self.visitRnemdblock(ctx.rnemdblock()))
        elif ctx.lightblock():
            result['light'].append(self.visitLightblock(ctx.lightblock()))
        elif ctx.minimizerblock():
            result['minimizer'].append(self.visitMinimizerblock(ctx.minimizerblock()))

    # -------------------------------------------------------------------------
    # Assignment
    # -------------------------------------------------------------------------

    def visitAssignment(self, ctx):
        """Returns (key, value) tuple."""
        key = ctx.ID().getText()
        value = self.visitConstant(ctx.constant())
        return key, value

    def visitConstant(self, ctx):
        if ctx.intConst():
            return self.visitIntConst(ctx.intConst())
        elif ctx.floatConst():
            return self.visitFloatConst(ctx.floatConst())
        elif ctx.vectorConst():
            return self.visitVectorConst(ctx.vectorConst())
        elif ctx.TRUE():
            return True
        elif ctx.FALSE():
            return False
        elif ctx.ID():
            return ctx.ID().getText()
        elif ctx.StringLiteral():
            # Strip surrounding quotes
            return ctx.StringLiteral().getText()[1:-1]
        return None

    # -------------------------------------------------------------------------
    # Simple blocks (component, zconstraint, restraint, flucq, rnemd,
    #                light, minimizer) — all share the same structure
    # -------------------------------------------------------------------------

    def _visit_simple_block(self, ctx):
        """Handle blocks that contain only assignments."""
        result = {'assignments': {}}
        for a in ctx.assignment():
            key, value = self.visitAssignment(a)
            result['assignments'][key] = value
        return result

    def visitComponentblock(self, ctx):
        return self._visit_simple_block(ctx)

    def visitZconstraintblock(self, ctx):
        return self._visit_simple_block(ctx)

    def visitRestraintblock(self, ctx):
        return self._visit_simple_block(ctx)

    def visitFlucqblock(self, ctx):
        return self._visit_simple_block(ctx)

    def visitRnemdblock(self, ctx):
        return self._visit_simple_block(ctx)

    def visitLightblock(self, ctx):
        return self._visit_simple_block(ctx)

    def visitMinimizerblock(self, ctx):
        return self._visit_simple_block(ctx)

    # -------------------------------------------------------------------------
    # Molecule block
    # -------------------------------------------------------------------------

    def visitMoleculeblock(self, ctx):
        result = {
            'assignments':  {},
            'atoms':        [],
            'bonds':        [],
            'bends':        [],
            'torsions':     [],
            'inversions':   [],
            'rigidbodies':  [],
            'cutoffgroups': [],
            'constraints':  [],
            'sequences':    [],
        }
        for stmt in ctx.moleculestatement():
            self._dispatch_moleculestatement(stmt, result)
        return result

    def _dispatch_moleculestatement(self, ctx, result):
        if ctx.assignment():
            key, value = self.visitAssignment(ctx.assignment())
            result['assignments'][key] = value
        elif ctx.atomblock():
            result['atoms'].append(self.visitAtomblock(ctx.atomblock()))
        elif ctx.bondblock():
            result['bonds'].append(self.visitBondblock(ctx.bondblock()))
        elif ctx.bendblock():
            result['bends'].append(self.visitBendblock(ctx.bendblock()))
        elif ctx.torsionblock():
            result['torsions'].append(self.visitTorsionblock(ctx.torsionblock()))
        elif ctx.inversionblock():
            result['inversions'].append(self.visitInversionblock(ctx.inversionblock()))
        elif ctx.rigidbodyblock():
            result['rigidbodies'].append(self.visitRigidbodyblock(ctx.rigidbodyblock()))
        elif ctx.cutoffgroupblock():
            result['cutoffgroups'].append(self.visitCutoffgroupblock(ctx.cutoffgroupblock()))
        elif ctx.constraintblock():
            result['constraints'].append(self.visitConstraintblock(ctx.constraintblock()))
        elif ctx.sequencestring():
            result['sequences'].append(self.visitSequencestring(ctx.sequencestring()))

    # -------------------------------------------------------------------------
    # Atom block
    # -------------------------------------------------------------------------

    def visitAtomblock(self, ctx):
        result = {
            'index':        self.visitIntConst(ctx.intConst()),
            'assignments':  {},
            'positions':    [],
            'orientations': [],
            'charges':      [],
        }
        for stmt in ctx.atomstatement():
            self._dispatch_atomstatement(stmt, result)
        return result

    def _dispatch_atomstatement(self, ctx, result):
        if ctx.assignment():
            key, value = self.visitAssignment(ctx.assignment())
            result['assignments'][key] = value
        elif ctx.POSITION():
            result['positions'].append(self.visitDoubleNumberTuple(ctx.doubleNumberTuple()))
        elif ctx.ORIENTATION():
            result['orientations'].append(self.visitDoubleNumberTuple(ctx.doubleNumberTuple()))
        elif ctx.CHARGE():
            result['charges'].append(self.visitFloatConst(ctx.floatConst()))

    # -------------------------------------------------------------------------
    # Bond block
    # -------------------------------------------------------------------------

    def visitBondblock(self, ctx):
        result = {
            'index':       self.visitIntConst(ctx.intConst()) if ctx.intConst() else None,
            'assignments': {},
            'members':     [],
            'potentials':  [],
        }
        for stmt in ctx.bondstatement():
            self._dispatch_bondstatement(stmt, result)
        return result

    def _dispatch_bondstatement(self, ctx, result):
        if ctx.assignment():
            key, value = self.visitAssignment(ctx.assignment())
            result['assignments'][key] = value
        elif ctx.MEMBERS():
            result['members'].append(self.visitInttuple(ctx.inttuple()))
        else:
            # Potential type — find which keyword token is present
            pot_type = self._bond_potential_type(ctx)
            if pot_type:
                result['potentials'].append({
                    'type':   pot_type,
                    'params': self.visitDoubleNumberTuple(ctx.doubleNumberTuple())
                              if ctx.doubleNumberTuple() else
                              [self.visitFloatConst(ctx.floatConst())]
                })

    def _bond_potential_type(self, ctx):
        for name in ('FIXED', 'HARMONIC', 'CUBIC', 'QUARTIC', 'POLYNOMIAL', 'MORSE'):
            if getattr(ctx, name, lambda: None)():
                return name
        return None

    # -------------------------------------------------------------------------
    # Bend block
    # -------------------------------------------------------------------------

    def visitBendblock(self, ctx):
        result = {
            'index':       self.visitIntConst(ctx.intConst()) if ctx.intConst() else None,
            'assignments': {},
            'members':     [],
            'potentials':  [],
        }
        for stmt in ctx.bendstatement():
            self._dispatch_bendstatement(stmt, result)
        return result

    def _dispatch_bendstatement(self, ctx, result):
        if ctx.assignment():
            key, value = self.visitAssignment(ctx.assignment())
            result['assignments'][key] = value
        elif ctx.MEMBERS():
            result['members'].append(self.visitInttuple(ctx.inttuple()))
        else:
            pot_type = self._bend_potential_type(ctx)
            if pot_type:
                result['potentials'].append({
                    'type':   pot_type,
                    'params': self.visitDoubleNumberTuple(ctx.doubleNumberTuple())
                })

    def _bend_potential_type(self, ctx):
        for name in ('HARMONIC', 'GHOSTBEND', 'UREYBRADLEY', 'CUBIC',
                     'QUARTIC', 'POLYNOMIAL', 'COSINE'):
            if getattr(ctx, name, lambda: None)():
                return name
        return None

    # -------------------------------------------------------------------------
    # Torsion block
    # -------------------------------------------------------------------------

    def visitTorsionblock(self, ctx):
        result = {
            'index':       self.visitIntConst(ctx.intConst()) if ctx.intConst() else None,
            'assignments': {},
            'members':     [],
            'potentials':  [],
        }
        for stmt in ctx.torsionstatement():
            self._dispatch_torsionstatement(stmt, result)
        return result

    def _dispatch_torsionstatement(self, ctx, result):
        if ctx.assignment():
            key, value = self.visitAssignment(ctx.assignment())
            result['assignments'][key] = value
        elif ctx.MEMBERS():
            result['members'].append(self.visitInttuple(ctx.inttuple()))
        else:
            pot_type = self._torsion_potential_type(ctx)
            if pot_type:
                result['potentials'].append({
                    'type':   pot_type,
                    'params': self.visitDoubleNumberTuple(ctx.doubleNumberTuple())
                })

    def _torsion_potential_type(self, ctx):
        for name in ('GHOSTTORSION', 'CUBIC', 'QUARTIC', 'POLYNOMIAL',
                     'CHARMM', 'OPLS', 'TRAPPE', 'HARMONIC'):
            if getattr(ctx, name, lambda: None)():
                return name
        return None

    # -------------------------------------------------------------------------
    # Inversion block
    # -------------------------------------------------------------------------

    def visitInversionblock(self, ctx):
        result = {
            'index':       self.visitIntConst(ctx.intConst()) if ctx.intConst() else None,
            'assignments': {},
            'center':      None,
            'satellites':  [],
            'potentials':  [],
        }
        for stmt in ctx.inversionstatement():
            self._dispatch_inversionstatement(stmt, result)
        return result

    def _dispatch_inversionstatement(self, ctx, result):
        if ctx.assignment():
            key, value = self.visitAssignment(ctx.assignment())
            result['assignments'][key] = value
        elif ctx.CENTER():
            result['center'] = self.visitIntConst(ctx.intConst())
        elif ctx.SATELLITES():
            result['satellites'] = self.visitInttuple(ctx.inttuple())
        else:
            pot_type = self._inversion_potential_type(ctx)
            if pot_type:
                result['potentials'].append({
                    'type':   pot_type,
                    'params': self.visitDoubleNumberTuple(ctx.doubleNumberTuple())
                })

    def _inversion_potential_type(self, ctx):
        for name in ('AMBERIMPROPER', 'IMPROPERCOSINE', 'HARMONIC',
                     'CENTRALATOMHEIGHT', 'DREIDING'):
            if getattr(ctx, name, lambda: None)():
                return name
        return None

    # -------------------------------------------------------------------------
    # RigidBody block
    # -------------------------------------------------------------------------

    def visitRigidbodyblock(self, ctx):
        result = {
            'index':       self.visitIntConst(ctx.intConst()),
            'assignments': {},
            'members':     [],
        }
        for stmt in ctx.rigidbodystatement():
            if stmt.assignment():
                key, value = self.visitAssignment(stmt.assignment())
                result['assignments'][key] = value
            elif stmt.MEMBERS():
                result['members'].append(self.visitInttuple(stmt.inttuple()))
        return result

    # -------------------------------------------------------------------------
    # CutoffGroup block
    # -------------------------------------------------------------------------

    def visitCutoffgroupblock(self, ctx):
        result = {
            'index':       self.visitIntConst(ctx.intConst()) if ctx.intConst() else None,
            'assignments': {},
            'members':     [],
        }
        for stmt in ctx.cutoffgroupstatement():
            if stmt.assignment():
                key, value = self.visitAssignment(stmt.assignment())
                result['assignments'][key] = value
            elif stmt.MEMBERS():
                result['members'].append(self.visitInttuple(stmt.inttuple()))
        return result

    # -------------------------------------------------------------------------
    # Constraint block
    # -------------------------------------------------------------------------

    def visitConstraintblock(self, ctx):
        result = {
            'index':       self.visitIntConst(ctx.intConst()) if ctx.intConst() else None,
            'assignments': {},
            'members':     [],
        }
        for stmt in ctx.constraintstatement():
            if stmt.assignment():
                key, value = self.visitAssignment(stmt.assignment())
                result['assignments'][key] = value
            elif stmt.MEMBERS():
                result['members'].append(self.visitInttuple(stmt.inttuple()))
        return result

    # -------------------------------------------------------------------------
    # Nodes block
    # -------------------------------------------------------------------------

    def visitNodesblock(self, ctx):
        result = {
            'index':       self.visitIntConst(ctx.intConst()) if ctx.intConst() else None,
            'assignments': {},
            'members':     [],
        }
        for stmt in ctx.nodesstatement():
            if stmt.assignment():
                key, value = self.visitAssignment(stmt.assignment())
                result['assignments'][key] = value
            elif stmt.MEMBERS():
                result['members'].append(self.visitInttuple(stmt.inttuple()))
        return result

    # -------------------------------------------------------------------------
    # Fragment block
    # -------------------------------------------------------------------------

    def visitFragmentblock(self, ctx):
        result = {
            'assignments':  {},
            'atoms':        [],
            'bonds':        [],
            'bends':        [],
            'torsions':     [],
            'inversions':   [],
            'rigidbodies':  [],
            'cutoffgroups': [],
            'constraints':  [],
            'nodes':        [],
        }
        for stmt in ctx.fragmentstatement():
            if stmt.assignment():
                key, value = self.visitAssignment(stmt.assignment())
                result['assignments'][key] = value
            elif stmt.atomblock():
                result['atoms'].append(self.visitAtomblock(stmt.atomblock()))
            elif stmt.bondblock():
                result['bonds'].append(self.visitBondblock(stmt.bondblock()))
            elif stmt.bendblock():
                result['bends'].append(self.visitBendblock(stmt.bendblock()))
            elif stmt.torsionblock():
                result['torsions'].append(self.visitTorsionblock(stmt.torsionblock()))
            elif stmt.inversionblock():
                result['inversions'].append(self.visitInversionblock(stmt.inversionblock()))
            elif stmt.rigidbodyblock():
                result['rigidbodies'].append(self.visitRigidbodyblock(stmt.rigidbodyblock()))
            elif stmt.cutoffgroupblock():
                result['cutoffgroups'].append(self.visitCutoffgroupblock(stmt.cutoffgroupblock()))
            elif stmt.constraintblock():
                result['constraints'].append(self.visitConstraintblock(stmt.constraintblock()))
            elif stmt.nodesblock():
                result['nodes'].append(self.visitNodesblock(stmt.nodesblock()))
        return result

    # -------------------------------------------------------------------------
    # Sequence string
    # -------------------------------------------------------------------------

    def visitSequencestring(self, ctx):
        return self.visitConstant(ctx.constant())

    # -------------------------------------------------------------------------
    # Numeric helpers
    # -------------------------------------------------------------------------

    def visitIntConst(self, ctx):
        text = ctx.getText()
        # Strip trailing L/l suffix
        stripped = text.rstrip('lL')
        try:
            return int(stripped, 0)  # base 0 handles hex/octal/decimal
        except ValueError:
            return text

    def visitFloatConst(self, ctx):
        return self._parse_float(ctx.getText())

    def visitDoubleNumber(self, ctx):
        if ctx.intConst():
            return self.visitIntConst(ctx.intConst())
        return self.visitFloatConst(ctx.floatConst())

    def visitDoubleNumberTuple(self, ctx):
        return [self.visitDoubleNumber(dn) for dn in ctx.doubleNumber()]

    def visitInttuple(self, ctx):
        return [self.visitIntConst(ic) for ic in ctx.intConst()]

    def visitVectorConst(self, ctx):
        return self.visitDoubleNumberTuple(ctx.doubleNumberTuple())

    def _parse_float(self, text):
        """
        Parse a float/double literal, handling Fortran-style exponent (d/D).
        e.g. 1.2d-04  ->  1.2e-04
        """
        # Strip trailing type suffix f/F/d/D that isn't part of an exponent
        # Replace Fortran exponent marker d/D with e for Python float()
        normalized = text.rstrip('fF')
        # Replace a leading d/D exponent marker (after digits) with e
        # e.g. 1.2d-4 -> 1.2e-4, but not a trailing bare 'd' suffix
        import re
        normalized = re.sub(r'([0-9])([dD])([+\-]?[0-9])', r'\1e\3', normalized)
        # Strip any remaining trailing d/D suffix (bare type marker)
        normalized = normalized.rstrip('dD')
        try:
            return float(normalized)
        except ValueError:
            return text
