# Generated from OMD.g4 by ANTLR 4.13.2
from antlr4 import *
if "." in __name__:
    from .OMDParser import OMDParser
else:
    from OMDParser import OMDParser

# This class defines a complete generic visitor for a parse tree produced by OMDParser.

class OMDVisitor(ParseTreeVisitor):

    # Visit a parse tree produced by OMDParser#omdfile.
    def visitOmdfile(self, ctx:OMDParser.OmdfileContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#statement.
    def visitStatement(self, ctx:OMDParser.StatementContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#assignment.
    def visitAssignment(self, ctx:OMDParser.AssignmentContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#constant.
    def visitConstant(self, ctx:OMDParser.ConstantContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#componentblock.
    def visitComponentblock(self, ctx:OMDParser.ComponentblockContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#zconstraintblock.
    def visitZconstraintblock(self, ctx:OMDParser.ZconstraintblockContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#restraintblock.
    def visitRestraintblock(self, ctx:OMDParser.RestraintblockContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#flucqblock.
    def visitFlucqblock(self, ctx:OMDParser.FlucqblockContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#rnemdblock.
    def visitRnemdblock(self, ctx:OMDParser.RnemdblockContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#lightblock.
    def visitLightblock(self, ctx:OMDParser.LightblockContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#minimizerblock.
    def visitMinimizerblock(self, ctx:OMDParser.MinimizerblockContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#moleculeblock.
    def visitMoleculeblock(self, ctx:OMDParser.MoleculeblockContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#moleculestatement.
    def visitMoleculestatement(self, ctx:OMDParser.MoleculestatementContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#atomblock.
    def visitAtomblock(self, ctx:OMDParser.AtomblockContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#atomstatement.
    def visitAtomstatement(self, ctx:OMDParser.AtomstatementContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#bondblock.
    def visitBondblock(self, ctx:OMDParser.BondblockContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#bondstatement.
    def visitBondstatement(self, ctx:OMDParser.BondstatementContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#bendblock.
    def visitBendblock(self, ctx:OMDParser.BendblockContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#bendstatement.
    def visitBendstatement(self, ctx:OMDParser.BendstatementContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#torsionblock.
    def visitTorsionblock(self, ctx:OMDParser.TorsionblockContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#torsionstatement.
    def visitTorsionstatement(self, ctx:OMDParser.TorsionstatementContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#inversionblock.
    def visitInversionblock(self, ctx:OMDParser.InversionblockContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#inversionstatement.
    def visitInversionstatement(self, ctx:OMDParser.InversionstatementContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#rigidbodyblock.
    def visitRigidbodyblock(self, ctx:OMDParser.RigidbodyblockContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#rigidbodystatement.
    def visitRigidbodystatement(self, ctx:OMDParser.RigidbodystatementContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#cutoffgroupblock.
    def visitCutoffgroupblock(self, ctx:OMDParser.CutoffgroupblockContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#cutoffgroupstatement.
    def visitCutoffgroupstatement(self, ctx:OMDParser.CutoffgroupstatementContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#nodesblock.
    def visitNodesblock(self, ctx:OMDParser.NodesblockContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#nodesstatement.
    def visitNodesstatement(self, ctx:OMDParser.NodesstatementContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#fragmentblock.
    def visitFragmentblock(self, ctx:OMDParser.FragmentblockContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#fragmentstatement.
    def visitFragmentstatement(self, ctx:OMDParser.FragmentstatementContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#constraintblock.
    def visitConstraintblock(self, ctx:OMDParser.ConstraintblockContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#constraintstatement.
    def visitConstraintstatement(self, ctx:OMDParser.ConstraintstatementContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#sequencestring.
    def visitSequencestring(self, ctx:OMDParser.SequencestringContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#doubleNumberTuple.
    def visitDoubleNumberTuple(self, ctx:OMDParser.DoubleNumberTupleContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#inttuple.
    def visitInttuple(self, ctx:OMDParser.InttupleContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#intConst.
    def visitIntConst(self, ctx:OMDParser.IntConstContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#doubleNumber.
    def visitDoubleNumber(self, ctx:OMDParser.DoubleNumberContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#floatConst.
    def visitFloatConst(self, ctx:OMDParser.FloatConstContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by OMDParser#vectorConst.
    def visitVectorConst(self, ctx:OMDParser.VectorConstContext):
        return self.visitChildren(ctx)



del OMDParser