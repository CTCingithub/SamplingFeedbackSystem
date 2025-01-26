import sympy as p
import itertools


def Translation_4x4(Displacement):
    # Calculate Transformation Matrix in translation cases
    return (
        p.Identity(3)
        .as_explicit()
        .row_join(Displacement)
        .col_join(p.Matrix([[0, 0, 0, 1]]))
    )


def Rotation_RPY_4x4(RPY):
    # psi, theta, phi = RPY[0], RPY[1], RPY[2]
    return (
        MatrixExp_4x4(p.Matrix([1, 0, 0, 0, 0, 0]), RPY[0])
        @ MatrixExp_4x4(p.Matrix([0, 1, 0, 0, 0, 0]), RPY[1])
        @ MatrixExp_4x4(p.Matrix([0, 0, 1, 0, 0, 0]), RPY[2])
    )


def TransformationMatrix_Inverse(TransformationMatrix):
    R_matrix = TransformationMatrix[:3, :3]
    p_vector = TransformationMatrix[:3, 3].reshape(3, 1)
    return R_matrix.T.row_join(-R_matrix.T @ p_vector).col_join(
        p.Matrix([[0, 0, 0, 1]])
    )


def Vector2Matrix_3x3(Vector):
    # Convert 3x1-shaped \hat{\omega} vectors to 3x3-shaped
    # [\hat{\omega}] matrixes
    #! Sympy doesn't support reshapes like .reshape(-1,?)
    vec = Vector.reshape(3, 1)
    return p.Matrix(
        [
            [0, -vec[2, 0], vec[1, 0]],
            [vec[2, 0], 0, -vec[0, 0]],
            [-vec[1, 0], vec[0, 0], 0],
        ]
    )


def MatrixExp_3x3(Vector, Angle):
    # Calculate exp([\hat{\omega}] \theta) matrixes, \hat{\omega} is a
    # 3x1-shaped vector
    mat_temp = Vector2Matrix_3x3(Vector)
    return (
        p.Identity(3)
        + mat_temp * p.sin(Angle)
        + mat_temp @ mat_temp * (1 - p.cos(Angle))
    ).as_explicit()


def Joint2Twist(Joint, Location):
    # Calculate twists from joint information
    return Joint.col_join(Vector2Matrix_3x3(Location) @ Joint)


def PVector(Twist, Angle):
    # Calculate \vec{p}s
    omega = p.Matrix(Twist[:3, :]).reshape(3, 1)
    v = p.Matrix(Twist[3:, :]).reshape(3, 1)
    return (p.Identity(3).as_explicit() - MatrixExp_3x3(omega, Angle)) @ (
        Vector2Matrix_3x3(omega) @ v
    ) + omega @ omega.T @ v * Angle


def MatrixExp_4x4(Twist, Angle):
    # Calculate exp([\hat{\xi}] \theta) matrixes, \hat{\xi} is a
    # 6x1-shaped twist
    UpperLeft = MatrixExp_3x3(Twist[:3, :], Angle)
    UpperRight = PVector(Twist, Angle)
    return UpperLeft.row_join(UpperRight).col_join(p.Matrix([[0, 0, 0, 1]]))


def Twist2Matrix_4x4(Twist):
    return (
        Vector2Matrix_3x3(Twist[:3, :])
        .row_join(Twist[3:, :])
        .col_join(p.Matrix([[0, 0, 0, 0]]))
    )


def AdjointMatrix(TransMatrix):
    R = TransMatrix[:3, :3]
    P = TransMatrix[:-1, -1]
    return p.simplify(
        R.row_join(p.zeros(3, 3)).col_join((Vector2Matrix_3x3(P) @ R).row_join(R))
    )


def AdjointInverseMatrix(TransMatrix):
    R = TransMatrix[:3, :3]
    P = TransMatrix[:-1, -1]
    return p.simplify(
        R.T.row_join(p.zeros(3, 3)).col_join(
            (-R.T @ Vector2Matrix_3x3(P)).row_join(R.T)
        )
    )


def MatlabCode(Matrix):
    # Generate Matlab code of matrixes
    Code = ""
    Mat = p.simplify(Matrix)
    for i in range(Mat.shape[0]):
        for j in range(Mat.shape[1]):
            if j == 0:
                Code = f"{Code}     {p.octave_code(Mat[i, j])}"
            else:
                Code = f"{Code}, {p.octave_code(Mat[i, j])}"
        Code = f"{Code}\n    " if i == Mat.shape[0] - 1 else f"{Code};\n   "
    return f"[\n    {Code[1:-2]}      ]"


def KMatrix_3x6(Omega):
    omega_x, omega_y, omega_z = Omega[0, 0], Omega[1, 0], Omega[2, 0]
    return p.Matrix(
        [
            [omega_x, omega_y, omega_z, 0, 0, 0],
            [0, omega_x, 0, omega_y, omega_z, 0],
            [0, 0, omega_x, 0, omega_y, omega_z],
        ]
    )


def Ttilde_10x1(SpatialVelocity, TransformationMat):
    return (
        (SpatialVelocity[:3, :].T / 2 @ KMatrix_3x6(SpatialVelocity[:3, :]))
        .row_join(
            SpatialVelocity[3:, :].T
            @ Vector2Matrix_3x3(SpatialVelocity[:3, :])
            @ TransformationMat[:3, :3]
        )
        .row_join(SpatialVelocity[3:, :].T / 2 @ SpatialVelocity[3:, :])
    ).T


def Vtilde_10x1(rVec, gVec, TransformationMat):
    return (
        p.zeros(1, 6)
        .row_join(-gVec.T @ TransformationMat[:3, :3])
        .row_join(-gVec.T @ rVec)
    ).T


def VecDiff(Vec1, Vec2):
    VecA = Vec1.reshape(len(Vec1), 1)
    VecB = Vec2.reshape(len(Vec2), 1)
    return VecA.jacobian(VecB).T


def VariablesGen(Variables):
    Code = "".join(f"{p.octave_code(i)}, " for i in Variables)
    return Code[:-2]


def CreateMatlabFunction(FunNameWithPath, FunName, Fun, Variables):
    FileName = f"{FunNameWithPath}.m"
    with open(FileName, "w") as File:
        File.write(f"function {FunName} = {FunName}(States)\n")
        File.write("    StateCell = num2cell(States);\n")
        File.write(f"    [{VariablesGen(Variables)}]" + " = deal(StateCell{:});\n")
        File.write(f"    {FunName} = ")
        File.write(MatlabCode(Fun))
        File.write(";\nend\n")


def RMat2YMat(R, v_vec, a_vec, phi_vec, psi_vec):
    # * Reshape state vectors
    V_Vec = v_vec.reshape(R.shape[0], 1)
    A_Vec = a_vec.reshape(R.shape[0], 1)
    Phi_Vec = phi_vec.reshape(R.shape[0], 1)
    Psi_Vec = psi_vec.reshape(R.shape[0], 1)

    # * Initialize Y matrix
    Y = p.zeros(R.shape[0], R.shape[1])
    for i, j in itertools.product(range(R.shape[0]), range(R.shape[1])):
        # ? Row index is i, Column index is j
        Y[i, j] = R[i, j]
        # ! Substitue acceleration relevant terms
        for acceleration_index in range(R.shape[0]):
            Y[i, j] = Y[i, j].subs(
                A_Vec[acceleration_index], Psi_Vec[acceleration_index]
            )
        # ! Substitue velocity relevant terms
        for velocity_index_1, velocity_index_2 in itertools.product(
            range(R.shape[0]), range(R.shape[0])
        ):
            Y[i, j] = Y[i, j].subs(
                V_Vec[velocity_index_1] * V_Vec[velocity_index_2],
                0.5
                * (
                    V_Vec[velocity_index_1] * Phi_Vec[velocity_index_2]
                    + V_Vec[velocity_index_2] * Phi_Vec[velocity_index_1]
                ),
            )

    return Y
