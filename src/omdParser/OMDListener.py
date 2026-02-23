# Generated from OMD.g4 by ANTLR 4.13.2
from antlr4 import *
if "." in __name__:
    from .OMDParser import OMDParser
else:
    from OMDParser import OMDParser

# This class defines a complete listener for a parse tree produced by OMDParser.
class OMDListener(ParseTreeListener):

    # Enter a parse tree produced by OMDParser#omdfile.
    def enterOmdfile(self, ctx:OMDParser.OmdfileContext):
        pass

    # Exit a parse tree produced by OMDParser#omdfile.
    def exitOmdfile(self, ctx:OMDParser.OmdfileContext):
        pass


    # Enter a parse tree produced by OMDParser#statement.
    def enterStatement(self, ctx:OMDParser.StatementContext):
        pass

    # Exit a parse tree produced by OMDParser#statement.
    def exitStatement(self, ctx:OMDParser.StatementContext):
        pass


    # Enter a parse tree produced by OMDParser#assignment.
    def enterAssignment(self, ctx:OMDParser.AssignmentContext):
        pass

    # Exit a parse tree produced by OMDParser#assignment.
    def exitAssignment(self, ctx:OMDParser.AssignmentContext):
        pass


    # Enter a parse tree produced by OMDParser#constant.
    def enterConstant(self, ctx:OMDParser.ConstantContext):
        pass

    # Exit a parse tree produced by OMDParser#constant.
    def exitConstant(self, ctx:OMDParser.ConstantContext):
        pass


    # Enter a parse tree produced by OMDParser#componentblock.
    def enterComponentblock(self, ctx:OMDParser.ComponentblockContext):
        pass

    # Exit a parse tree produced by OMDParser#componentblock.
    def exitComponentblock(self, ctx:OMDParser.ComponentblockContext):
        pass


    # Enter a parse tree produced by OMDParser#zconstraintblock.
    def enterZconstraintblock(self, ctx:OMDParser.ZconstraintblockContext):
        pass

    # Exit a parse tree produced by OMDParser#zconstraintblock.
    def exitZconstraintblock(self, ctx:OMDParser.ZconstraintblockContext):
        pass


    # Enter a parse tree produced by OMDParser#restraintblock.
    def enterRestraintblock(self, ctx:OMDParser.RestraintblockContext):
        pass

    # Exit a parse tree produced by OMDParser#restraintblock.
    def exitRestraintblock(self, ctx:OMDParser.RestraintblockContext):
        pass


    # Enter a parse tree produced by OMDParser#flucqblock.
    def enterFlucqblock(self, ctx:OMDParser.FlucqblockContext):
        pass

    # Exit a parse tree produced by OMDParser#flucqblock.
    def exitFlucqblock(self, ctx:OMDParser.FlucqblockContext):
        pass


    # Enter a parse tree produced by OMDParser#rnemdblock.
    def enterRnemdblock(self, ctx:OMDParser.RnemdblockContext):
        pass

    # Exit a parse tree produced by OMDParser#rnemdblock.
    def exitRnemdblock(self, ctx:OMDParser.RnemdblockContext):
        pass


    # Enter a parse tree produced by OMDParser#lightblock.
    def enterLightblock(self, ctx:OMDParser.LightblockContext):
        pass

    # Exit a parse tree produced by OMDParser#lightblock.
    def exitLightblock(self, ctx:OMDParser.LightblockContext):
        pass


    # Enter a parse tree produced by OMDParser#minimizerblock.
    def enterMinimizerblock(self, ctx:OMDParser.MinimizerblockContext):
        pass

    # Exit a parse tree produced by OMDParser#minimizerblock.
    def exitMinimizerblock(self, ctx:OMDParser.MinimizerblockContext):
        pass


    # Enter a parse tree produced by OMDParser#moleculeblock.
    def enterMoleculeblock(self, ctx:OMDParser.MoleculeblockContext):
        pass

    # Exit a parse tree produced by OMDParser#moleculeblock.
    def exitMoleculeblock(self, ctx:OMDParser.MoleculeblockContext):
        pass


    # Enter a parse tree produced by OMDParser#moleculestatement.
    def enterMoleculestatement(self, ctx:OMDParser.MoleculestatementContext):
        pass

    # Exit a parse tree produced by OMDParser#moleculestatement.
    def exitMoleculestatement(self, ctx:OMDParser.MoleculestatementContext):
        pass


    # Enter a parse tree produced by OMDParser#atomblock.
    def enterAtomblock(self, ctx:OMDParser.AtomblockContext):
        pass

    # Exit a parse tree produced by OMDParser#atomblock.
    def exitAtomblock(self, ctx:OMDParser.AtomblockContext):
        pass


    # Enter a parse tree produced by OMDParser#atomstatement.
    def enterAtomstatement(self, ctx:OMDParser.AtomstatementContext):
        pass

    # Exit a parse tree produced by OMDParser#atomstatement.
    def exitAtomstatement(self, ctx:OMDParser.AtomstatementContext):
        pass


    # Enter a parse tree produced by OMDParser#bondblock.
    def enterBondblock(self, ctx:OMDParser.BondblockContext):
        pass

    # Exit a parse tree produced by OMDParser#bondblock.
    def exitBondblock(self, ctx:OMDParser.BondblockContext):
        pass


    # Enter a parse tree produced by OMDParser#bondstatement.
    def enterBondstatement(self, ctx:OMDParser.BondstatementContext):
        pass

    # Exit a parse tree produced by OMDParser#bondstatement.
    def exitBondstatement(self, ctx:OMDParser.BondstatementContext):
        pass


    # Enter a parse tree produced by OMDParser#bendblock.
    def enterBendblock(self, ctx:OMDParser.BendblockContext):
        pass

    # Exit a parse tree produced by OMDParser#bendblock.
    def exitBendblock(self, ctx:OMDParser.BendblockContext):
        pass


    # Enter a parse tree produced by OMDParser#bendstatement.
    def enterBendstatement(self, ctx:OMDParser.BendstatementContext):
        pass

    # Exit a parse tree produced by OMDParser#bendstatement.
    def exitBendstatement(self, ctx:OMDParser.BendstatementContext):
        pass


    # Enter a parse tree produced by OMDParser#torsionblock.
    def enterTorsionblock(self, ctx:OMDParser.TorsionblockContext):
        pass

    # Exit a parse tree produced by OMDParser#torsionblock.
    def exitTorsionblock(self, ctx:OMDParser.TorsionblockContext):
        pass


    # Enter a parse tree produced by OMDParser#torsionstatement.
    def enterTorsionstatement(self, ctx:OMDParser.TorsionstatementContext):
        pass

    # Exit a parse tree produced by OMDParser#torsionstatement.
    def exitTorsionstatement(self, ctx:OMDParser.TorsionstatementContext):
        pass


    # Enter a parse tree produced by OMDParser#inversionblock.
    def enterInversionblock(self, ctx:OMDParser.InversionblockContext):
        pass

    # Exit a parse tree produced by OMDParser#inversionblock.
    def exitInversionblock(self, ctx:OMDParser.InversionblockContext):
        pass


    # Enter a parse tree produced by OMDParser#inversionstatement.
    def enterInversionstatement(self, ctx:OMDParser.InversionstatementContext):
        pass

    # Exit a parse tree produced by OMDParser#inversionstatement.
    def exitInversionstatement(self, ctx:OMDParser.InversionstatementContext):
        pass


    # Enter a parse tree produced by OMDParser#rigidbodyblock.
    def enterRigidbodyblock(self, ctx:OMDParser.RigidbodyblockContext):
        pass

    # Exit a parse tree produced by OMDParser#rigidbodyblock.
    def exitRigidbodyblock(self, ctx:OMDParser.RigidbodyblockContext):
        pass


    # Enter a parse tree produced by OMDParser#rigidbodystatement.
    def enterRigidbodystatement(self, ctx:OMDParser.RigidbodystatementContext):
        pass

    # Exit a parse tree produced by OMDParser#rigidbodystatement.
    def exitRigidbodystatement(self, ctx:OMDParser.RigidbodystatementContext):
        pass


    # Enter a parse tree produced by OMDParser#cutoffgroupblock.
    def enterCutoffgroupblock(self, ctx:OMDParser.CutoffgroupblockContext):
        pass

    # Exit a parse tree produced by OMDParser#cutoffgroupblock.
    def exitCutoffgroupblock(self, ctx:OMDParser.CutoffgroupblockContext):
        pass


    # Enter a parse tree produced by OMDParser#cutoffgroupstatement.
    def enterCutoffgroupstatement(self, ctx:OMDParser.CutoffgroupstatementContext):
        pass

    # Exit a parse tree produced by OMDParser#cutoffgroupstatement.
    def exitCutoffgroupstatement(self, ctx:OMDParser.CutoffgroupstatementContext):
        pass


    # Enter a parse tree produced by OMDParser#nodesblock.
    def enterNodesblock(self, ctx:OMDParser.NodesblockContext):
        pass

    # Exit a parse tree produced by OMDParser#nodesblock.
    def exitNodesblock(self, ctx:OMDParser.NodesblockContext):
        pass


    # Enter a parse tree produced by OMDParser#nodesstatement.
    def enterNodesstatement(self, ctx:OMDParser.NodesstatementContext):
        pass

    # Exit a parse tree produced by OMDParser#nodesstatement.
    def exitNodesstatement(self, ctx:OMDParser.NodesstatementContext):
        pass


    # Enter a parse tree produced by OMDParser#fragmentblock.
    def enterFragmentblock(self, ctx:OMDParser.FragmentblockContext):
        pass

    # Exit a parse tree produced by OMDParser#fragmentblock.
    def exitFragmentblock(self, ctx:OMDParser.FragmentblockContext):
        pass


    # Enter a parse tree produced by OMDParser#fragmentstatement.
    def enterFragmentstatement(self, ctx:OMDParser.FragmentstatementContext):
        pass

    # Exit a parse tree produced by OMDParser#fragmentstatement.
    def exitFragmentstatement(self, ctx:OMDParser.FragmentstatementContext):
        pass


    # Enter a parse tree produced by OMDParser#constraintblock.
    def enterConstraintblock(self, ctx:OMDParser.ConstraintblockContext):
        pass

    # Exit a parse tree produced by OMDParser#constraintblock.
    def exitConstraintblock(self, ctx:OMDParser.ConstraintblockContext):
        pass


    # Enter a parse tree produced by OMDParser#constraintstatement.
    def enterConstraintstatement(self, ctx:OMDParser.ConstraintstatementContext):
        pass

    # Exit a parse tree produced by OMDParser#constraintstatement.
    def exitConstraintstatement(self, ctx:OMDParser.ConstraintstatementContext):
        pass


    # Enter a parse tree produced by OMDParser#sequencestring.
    def enterSequencestring(self, ctx:OMDParser.SequencestringContext):
        pass

    # Exit a parse tree produced by OMDParser#sequencestring.
    def exitSequencestring(self, ctx:OMDParser.SequencestringContext):
        pass


    # Enter a parse tree produced by OMDParser#doubleNumberTuple.
    def enterDoubleNumberTuple(self, ctx:OMDParser.DoubleNumberTupleContext):
        pass

    # Exit a parse tree produced by OMDParser#doubleNumberTuple.
    def exitDoubleNumberTuple(self, ctx:OMDParser.DoubleNumberTupleContext):
        pass


    # Enter a parse tree produced by OMDParser#inttuple.
    def enterInttuple(self, ctx:OMDParser.InttupleContext):
        pass

    # Exit a parse tree produced by OMDParser#inttuple.
    def exitInttuple(self, ctx:OMDParser.InttupleContext):
        pass


    # Enter a parse tree produced by OMDParser#intConst.
    def enterIntConst(self, ctx:OMDParser.IntConstContext):
        pass

    # Exit a parse tree produced by OMDParser#intConst.
    def exitIntConst(self, ctx:OMDParser.IntConstContext):
        pass


    # Enter a parse tree produced by OMDParser#doubleNumber.
    def enterDoubleNumber(self, ctx:OMDParser.DoubleNumberContext):
        pass

    # Exit a parse tree produced by OMDParser#doubleNumber.
    def exitDoubleNumber(self, ctx:OMDParser.DoubleNumberContext):
        pass


    # Enter a parse tree produced by OMDParser#floatConst.
    def enterFloatConst(self, ctx:OMDParser.FloatConstContext):
        pass

    # Exit a parse tree produced by OMDParser#floatConst.
    def exitFloatConst(self, ctx:OMDParser.FloatConstContext):
        pass


    # Enter a parse tree produced by OMDParser#vectorConst.
    def enterVectorConst(self, ctx:OMDParser.VectorConstContext):
        pass

    # Exit a parse tree produced by OMDParser#vectorConst.
    def exitVectorConst(self, ctx:OMDParser.VectorConstContext):
        pass



del OMDParser