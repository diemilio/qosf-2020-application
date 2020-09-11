# Import dependencies
import numpy as np
from qiskit import QuantumCircuit

def map_gates(qc_ori, trans_id = True, keep_H = False):

    """Define new circuit of same # of qubits as circuit passed (no classical bits, but can be added):"""
    qc_trans = QuantumCircuit(qc_ori.num_qubits)
    
    """Loop through all gates of original circuit and translate into basis gates: rx, rz, cz"""
    for i in range(len(qc_ori.data)):
        if type(qc_ori.data[i][0]).__name__ == 'RXGate':
            """RX is left unchanged. Rotation angle for RX gate is extracted from original circuit
               using .params method"""
            ind = qc_ori.data[i][1][0].index
            theta = qc_ori.data[i][0].params[0]
            qc_trans.rx(theta,ind)

        elif type(qc_ori.data[i][0]).__name__ == 'RZGate':
            """RZ is left unchanged. Rotation angle for RZ gate is extracted from original circuit
               using .params method"""
            ind = qc_ori.data[i][1][0].index
            phi = qc_ori.data[i][0].params[0]
            qc_trans.rz(phi,ind)
            
        elif type(qc_ori.data[i][0]).__name__ == 'RYGate':
            """RY is implemented using sequence of RX and RZ rotations as shown in the description above"""
            ind = qc_ori.data[i][1][0].index
            theta = qc_ori.data[i][0].params[0]
            qc_trans.rx(np.pi/2,ind)
            qc_trans.rz(theta,ind)
            qc_trans.rx(-np.pi/2,ind)
            
        elif type(qc_ori.data[i][0]).__name__ == 'IGate':
            """I is implemented as 0 rad z-rotation if trans_id is true. Otherwise, I is ignored"""
            if trans_id == True:
                ind = qc_ori.data[i][1][0].index
                theta = 0
                qc_trans.rz(theta,ind)

        elif type(qc_ori.data[i][0]).__name__ == 'XGate':
            """X is implemented by a pi rotation about the X axis"""
            ind = qc_ori.data[i][1][0].index
            theta = np.pi
            qc_trans.rx(theta,ind)
            
        elif type(qc_ori.data[i][0]).__name__ == 'ZGate':
            """Z is implemented by a pi rotation about the Z axis"""
            ind = qc_ori.data[i][1][0].index
            phi = np.pi
            qc_trans.rz(phi,ind)  

        elif type(qc_ori.data[i][0]).__name__ == 'YGate':
            """Y is implemented by a pi rotation about the Y axis"""
            ind = qc_ori.data[i][1][0].index
            theta = np.pi
            qc_trans.rx(np.pi/2,ind)
            qc_trans.rz(theta,ind)
            qc_trans.rx(-np.pi/2,ind)            
                  
        elif type(qc_ori.data[i][0]).__name__ == 'HGate':
            """H implemented as pi/2 rad z-rotation followed by  pi/2 rad x-rotation"""
            ind = qc_ori.data[i][1][0].index
            
            if keep_H:
                """If True, H gates are not translated into Rx, Rz basis"""
                qc_trans.h(ind)
            else:
                phi = np.pi
                theta = np.pi/2
                qc_trans.rz(phi,ind)
                qc_trans.rx(np.pi/2,ind)
                qc_trans.rz(theta,ind)
                qc_trans.rx(-np.pi/2,ind)

        elif type(qc_ori.data[i][0]).__name__ == 'CZGate':
            """CZ is left unchanged."""
            ind0 = qc_ori.data[i][1][0].index
            ind1 = qc_ori.data[i][1][1].index
            qc_trans.cz(ind0,ind1)            

        elif type(qc_ori.data[i][0]).__name__ == 'CXGate':
            """CX is replace by a CZ and Hadamard gates before and after CZ on the target qubit."""
            ind0 = qc_ori.data[i][1][0].index
            ind1 = qc_ori.data[i][1][1].index
            
            if keep_H:
                """If True, H gates are not translated into Rx, Rz basis"""
                qc_trans.h(ind1)
                qc_trans.cz(ind0,ind1)
                qc_trans.h(ind1)
            else:
                phi = np.pi
                theta = np.pi/2
                qc_trans.rz(phi,ind1)
                qc_trans.rx(np.pi/2,ind1)
                qc_trans.rz(theta,ind1)
                qc_trans.rx(-np.pi/2,ind1)
                qc_trans.cz(ind0,ind1)              
                qc_trans.rz(phi,ind1)
                qc_trans.rx(np.pi/2,ind1)
                qc_trans.rz(theta,ind1)
                qc_trans.rx(-np.pi/2,ind1)
            
        else:
            name = type(qc_ori.data[i][0]).__name__
            raise Exception('Circuit contains the unsupported gate: {}'.format(name))
            break
            
    return qc_trans

def optimize1(qc_ori):
    
    qc_red = qc_ori.copy()                                # Create copy of original circuit
    a = 1
    b = len(qc_red.data)
    reducible1 = ['XGate', 'YGate', 'ZGate', 'HGate']     # List of 1-qubit reducible gates
    twoqubitgt = ['CXGate', 'CZGate']                     # List of 2-qubit gates

    while a < b:
        gate_type1 = type(qc_red.data[a-1][0]).__name__   # Type of current gate

        """Check if gate is in list of 1-qubit reducible gates, and find if adjacent gate is of same type"""
        if gate_type1 in reducible1:
            gate_qb1 = qc_red.data[a-1][1][0].index             # qubit number of current gate

            mtch_flg = False
            j = a
            gate_qb2b = None

            while not mtch_flg:
                gate_type2 = type(qc_red.data[j][0]).__name__   # Type of gate being compared
                gate_qb2a = qc_red.data[j][1][0].index          # qubit number of gate being compared

                if gate_type2 in twoqubitgt:                    # Check qubit connections of 2-qubit gates
                    gate_qb2b = qc_red.data[j][1][1].index
                else:
                    gate_qb2b = None

                if (gate_qb2a == gate_qb1) or (gate_qb2b == gate_qb1):
                    if gate_type2 == gate_type1:
                        """Remove adjancent gates, get new length of circuit"""
                        qc_red.data.pop(j)
                        qc_red.data.pop(a-1)
                        b = len(qc_red.data)
                        a = 1
                    else:
                        a+=1

                    mtch_flg = True

                elif j < b-1:
                    j+=1
                else:
                    a+=1
                    mtch_flg = True

        elif gate_type1 == 'CZGate':
            """Check if gate is CZGate, and find if adjacent gate is of same type"""
            gate_qb1a = qc_red.data[a-1][1][0].index             # qubit number of control
            gate_qb1b = qc_red.data[a-1][1][1].index             # qubit number of target

            mtch_flg = False
            j = a
            gate_qb2b = None

            while not mtch_flg:
                gate_type2 = type(qc_red.data[j][0]).__name__   # Type of gate being compared
                gate_qb2a = qc_red.data[j][1][0].index          # first qubit gate being compared

                if gate_type2 in twoqubitgt:                    # Check qubit connections of 2-qubit gates
                    gate_qb2b = qc_red.data[j][1][1].index
                else:
                    gate_qb2b = None

                if (gate_qb2a == gate_qb1a) or (gate_qb2b == gate_qb1a):
                    if (gate_type2 == gate_type1 == 'CZGate') and (((gate_qb2a == gate_qb1b) and (gate_qb2b == gate_qb1a)) or ((gate_qb2a == gate_qb1b) and (gate_qb2b == gate_qb1a))):
                        """Check for adjacent CZ gates with matching or opposite qubits"""
                        qc_red.data.pop(j)
                        qc_red.data.pop(a-1)
                        b = len(qc_red.data)
                        a = 1
                    else:
                        a+=1

                    mtch_flg = True

                elif j < b-1:
                    j+=1
                else:
                    a+=1
                    mtch_flg = True
        else:
            a += 1
            
    return qc_red

def optimize2(qc_ori):
        
    qc_red = qc_ori.copy()                                # Create copy of original circuit
    a = 1
    b = len(qc_red.data)
    reducible = ['RXGate', 'RZGate']                      # List of 1-qubit reducible gates
    twoqubitgt = ['CXGate', 'CZGate']                     # List of 2-qubit gates

    while a < b:
        gate_type1 = type(qc_red.data[a-1][0]).__name__   # Type of current gate

        """Check if gate is in list of 1-qubit reducible gates, and find if adjacent gate is of same type"""
        if gate_type1 in reducible:
            gate_qb1 = qc_red.data[a-1][1][0].index             # qubit number of current gate
            gate_ang1 = qc_red.data[a-1][0].params[0]           # angle of rotation of current gate

            mtch_flg = False
            j = a
            gate_qb2b = None

            while not mtch_flg:
                gate_type2 = type(qc_red.data[j][0]).__name__   # Type of gate being compared
                gate_qb2a = qc_red.data[j][1][0].index          # qubit number of gate being compared

                if gate_type2 in twoqubitgt:                    # Check qubit connections of 2-qubit gates
                    gate_qb2b = qc_red.data[j][1][1].index
                else:
                    gate_qb2b = None

                if (gate_qb2a == gate_qb1) or (gate_qb2b == gate_qb1):
                    if gate_type2 == gate_type1:
                        """Add up rotation angles and remove one of the gates.
                           If sum of angles is 0, remove both gates"""
                        gate_ang2 = qc_red.data[j][0].params[0]         # angle of rotation of gate being compared
                        ang_sum = gate_ang2 + gate_ang1
                        qc_red.data[j][0].params[0] = ang_sum
                        
                        """If angle is multiple of 2*pi, remove gate"""
                        if ((ang_sum/(2*np.pi)) % 1) == 0:
                            qc_red.data.pop(j)
                        
                        qc_red.data.pop(a-1)
                        a = 1
                        b = len(qc_red.data)
                        
                    else:
                        a+=1

                    mtch_flg = True

                elif j < b-1:
                    j+=1
                else:
                    a+=1
                    mtch_flg = True
                    
        elif gate_type1 == 'CZGate':
            """Check if gate is CZGate, and find if adjacent gate is of same type"""
            gate_qb1a = qc_red.data[a-1][1][0].index             # qubit number of control
            gate_qb1b = qc_red.data[a-1][1][1].index             # qubit number of target

            mtch_flg = False
            j = a
            gate_qb2b = None

            while not mtch_flg:
                gate_type2 = type(qc_red.data[j][0]).__name__   # Type of gate being compared
                gate_qb2a = qc_red.data[j][1][0].index          # first qubit gate being compared

                if gate_type2 in twoqubitgt:                    # Check qubit connections of 2-qubit gates
                    gate_qb2b = qc_red.data[j][1][1].index
                else:
                    gate_qb2b = None

                if (gate_qb2a == gate_qb1a) or (gate_qb2b == gate_qb1a):
                    if (gate_type2 == gate_type1 == 'CZGate') and (((gate_qb2a == gate_qb1b) and (gate_qb2b == gate_qb1a)) or ((gate_qb2a == gate_qb1b) and (gate_qb2b == gate_qb1a))):
                        """Check for adjacent CZ gates with matching or opposite qubits"""
                        qc_red.data.pop(j)
                        qc_red.data.pop(a-1)
                        b = len(qc_red.data)
                        a = 1
                    else:
                        a+=1

                    mtch_flg = True

                elif j < b-1:
                    j+=1
                else:
                    a+=1
                    mtch_flg = True

        elif gate_type1 == 'HGate':
            """Check if gate is HGate, and find if adjacent gate is of same type"""
            gate_qb1 = qc_red.data[a-1][1][0].index             # qubit number of current gate

            mtch_flg = False
            j = a
            gate_qb2b = None

            while not mtch_flg:
                gate_type2 = type(qc_red.data[j][0]).__name__   # Type of gate being compared
                gate_qb2a = qc_red.data[j][1][0].index          # qubit number of gate being compared

                if gate_type2 in twoqubitgt:                    # Check qubit connections of 2-qubit gates
                    gate_qb2b = qc_red.data[j][1][1].index
                else:
                    gate_qb2b = None

                if (gate_qb2a == gate_qb1) or (gate_qb2b == gate_qb1):
                    if gate_type2 == gate_type1:
                        """Remove adjancent gates, get new length of circuit"""
                        qc_red.data.pop(j)
                        qc_red.data.pop(a-1)
                        b = len(qc_red.data)
                        a = 1
                    else:
                        a+=1

                    mtch_flg = True

                elif j < b-1:
                    j+=1
                else:
                    a+=1
                    mtch_flg = True            
            
        else:
            a += 1
            
    return qc_red

def translate(qc_ori, trans_id = False, op_level = 0):
    if op_level == 3:
        qc_trans3a = map_gates(qc_ori, trans_id, keep_H = True)
        qc_trans3b = optimize2(qc_trans3a)
        qc_trans3c = map_gates(qc_trans3b, trans_id, keep_H = False)
        qc_trans = optimize2(qc_trans3c)        
    
    elif op_level == 2:
        qc_trans2 = map_gates(qc_ori, trans_id)
        qc_trans = optimize2(qc_trans2)
    
    elif op_level == 1:
        qc_trans1 = optimize1(qc_ori)
        qc_trans = map_gates(qc_trans1, trans_id)
    
    elif op_level == 0:
        qc_trans = map_gates(qc_ori, trans_id)

    else:
        raise Exception('Optimization level must be 0, 1, 2 or 3')    
        
    return qc_trans