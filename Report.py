EPS = 1e-9

def matrixout(mx, size): #행렬(mx)와 행렬크기(size)를 출력
    for i in range(size):
        print("%d행" %(i+1)) #출력하는 행 구분
        for j in range(size):
            print("%7.3f" %mx[i][j], end=" ") #소수 3자리까지 각 행렬 성분을 출력
        print()
def transposeMatrix(m):
    return [[m[j][i] for j in range(len(m))] for i in range(len(m[0]))] #행렬(m)을 대각선을 기준으로 전치 후 반환

def getMatrixMinor(m, i ,j):
    return [row[:j] + row[j+1:] for row in (m[:i] + m[i+1:])] #i번째 행과 j번째 열을 제거한 소행렬을 반환

def getMatrixDeterminant(m):
    if len(m) == 1:
        return m[0][0] #행렬의 크기가 1인 경우
    if len(m) == 2:
        return m[0][0]*m[1][1] - m[0][1]*m[1][0] #행렬의 크기가 2인 경우

    determinant = 0
    for c in range(len(m)): #행렬의 크기가 3이상인 경우
        determinant += ((-1)**c)*m[0][c] * getMatrixDeterminant(getMatrixMinor(m,0,c)) # 재귀를 통해 det값을 구함
    return determinant


def getMatrixInverse(m): # det=0, 해당 함수를 호출하지 않음
    determinant = getMatrixDeterminant(m)

    if len(m) == 1:
        return [[1.0 / m[0][0]]] #행렬의 크기가 1인 경우
    if len(m) == 2:
        return [[m[1][1]/determinant, -1*m[0][1] / determinant],
                 [-1*m[1][0] / determinant, m[0][0] / determinant]] #행렬의 크기가 2인 경우

    cofactors = [] # 여인수 행렬
    for r in range(len(m)):
        cofactorRow = [] # 여인수 행렬의 행
        for c in range(len(m)):
            minor = getMatrixMinor(m,r,c) # 행렬의 소행렬을 대입
            cofactorRow.append(((-1)**(r+c)) * getMatrixDeterminant(minor)) # 여인수 행렬 행에 추가
        cofactors.append(cofactorRow) # 구한 여인수 행렬의 행을 여인수 행렬에 추가

    adjugate = transposeMatrix(cofactors) # 여인수 행렬을 전치 후 대입

    for r in range(len(adjugate)):
        for c in range(len(adjugate)):
            adjugate[r][c] = adjugate[r][c] / determinant # 수반행렬 / det = 역행렬을 이용

    return adjugate # 역행렬 반환

def has_inverse(m):
    determinant = getMatrixDeterminant(m) # det 대입
    if abs(determinant) < EPS: # 부동소수점 오차 고려하고 abs를 사용해서 det=0인 경우를 찾아냄
        print("행렬식이 0이므로 역행렬이 존재하지 않습니다. (행렬식)")
        return False
    return True # det =/= 0

#-------------------------------------가ㅏ우스-조던---------------------------------------------------------------------------------------------

def getMatrixInverse_by_gauss_jordan(m):
    n = len(m)
    aug = [list(m[i]) + [1.0 if i == j else 0.0 for j in range(n)] for i in range(n)] # 단위행렬 추가

    for col in range(n):
        pivot = col # 피벗은 해당 열에서 절댓값이 가장 큰 성분의 위치를 나타냄
        for r in range(col, n):
            if abs(aug[r][col]) > abs(aug[pivot][col]):
                pivot = r
        if abs(aug[pivot][col]) < EPS: # 절댓값이 가장 큰 피벗이 0인 경우를 찾아냄
            print("피벗이 0이므로 역행렬이 존재하지 않습니다. (가우스-조던)")
            return None

        if pivot != col: # col행을 기준으로 계산
            aug[col], aug[pivot] = aug[pivot], aug[col] # 피벗이 col행이 아닌 경우 위치 변환

        p = aug[col][col] # 피벗 실제 값(pivot은 피벗의 위치임)
        aug[col] = [x / p for x in aug[col]] # 피벗을 1로 만듦


        for r in range(n):
            if r != col and abs(aug[r][col]) > EPS: # 피벗행과 성분=0인 행 제외
                f = aug[r][col]
                new_row = []
                for rc in range(len(aug[r])): # 다른 행의 col 열을 0으로 만들기
                    a = aug[r][rc]             # 이 행의 값
                    b = aug[col][rc]           # 피벗 행의 같은 위치 값
                    new_row.append(a - f * b)
                aug[r] = new_row

    return [row[n:] for row in aug]

#-----------------------비교 및 검산-------------------------------------------
def compare_matrix(a, b):
    for i in range(len(a)):
        for j in range(len(a)):
            if abs(a[i][j] - b[i][j]) > EPS: # 부동소수점 오차 고려
                return False
    return True

def matrix_mul(a,b): #a=행렬, b=역행렬
    n=len(a)
    I = [[0.0] * n for t in range(n)] # 초기화
    for i in range(n):            
        for j in range(n):
            for k in range(n):
                I[i][j] += a[i][k]*b[k][j]; # 행렬과 역행렬의 곱셈
            if abs(I[i][j]) < EPS: # 부동소수점 오차 고려, 소수점 정리
                I[i][j] = 0.0
    return I

def identityMatrix(I1):
    n = len(I1)
    for i in range(n):
        for j in range(n):
            expected = 1.0 if i == j else 0.0  # 단위행렬은 대각 성분만 1임을 이용
            if abs(I1[i][j] - expected) > EPS: # 부동소수점 오차 고려
                return False
    return True
    
            

def main():
    while True:
        try:
            k = int(input("정방행렬의 차수 n을 입력하시오: ")) # 차수 입력
            if k <= 0:                                     # 차수가 0 이하인 경우 다시 반복
                print("차수는 양의 정수여야 합니다.")
                continue
            break
        except ValueError:                                  # 정수가 아닌 입력값을 받은 경우
            print("입력 오류: 정수를 입력해주세요.")

    print(f"{k}x{k} 행렬 A를 한 행씩 입력하세요 (값은 공백으로 구분):")
    matrix_a = []
    for i in range(k):
        while True:
            try:
                row_input = input(f"{i+1}행: ").strip()       # 문자열을 받고 공백 제거
                row_values = [float(x) for x in row_input.split()] # 문자열을 리스트로 변환
                if len(row_values) != k:                     # k개 입력하지 않은 경우 다시 반복
                    print(f"입력 오류: 정확히 {k}개의 값을 입력해야 합니다.")
                    continue
                matrix_a.append(row_values)                  # 입력받은 값을 행렬(리스트)에 추가
                break
            except ValueError:
                print("입력 오류: 숫자만 입력해주세요.")

    inverse=None    # 행렬식으로 구한 역행렬
    inverse2=None   # 가우스-조던방식으로 구한 역행렬

    #  행렬식으로 구한 역행렬
    try:
        if has_inverse(matrix_a): # det=/=0으로 역행렬 여부 판별
            print("\n행렬식으로 구한 역행렬:")
            inverse = getMatrixInverse(matrix_a)
            matrixout(inverse, k) # 역행렬 출력
    except Exception as e:
        pass

    # 가우스-조던 소거법으로 구한 역행렬
    try:
        inverse2 = getMatrixInverse_by_gauss_jordan(matrix_a) # 피벗으로 역행렬 여부 판별
        if inverse2 is not None:
            print("\n가우스-조던 소거법으로 구한 역행렬:")
            matrixout(inverse2, k) # 역행렬 출력
    except Exception as e:
        pass

    try:
        if inverse == None and inverse2 == None: # 역행렬이 존재하지 않으면 None
            print("\n두 방법 모두 역행렬이 존재하지 않습니다.")
        elif compare_matrix(inverse, inverse2): # 역행렬이 일치하는 경우
            print("\n두 방법의 결과가 동일합니다.")

            I1=matrix_mul(matrix_a,inverse)
            I2=matrix_mul(matrix_a,inverse2)

            if compare_matrix(I1, I2) and identityMatrix(I1) : # 역행렬이 맞는지 검산 후 검산 결과 출력
                print("\n검산 -> 행렬 x 역행렬 = 단위행렬")
                matrixout(I1, k)
            else:
                print("\n검산 결과 단위행렬이 아닙니다.")
        elif inverse is None or inverse2 is None:
            print("\n두 방법의 판정이 서로 다릅니다.")
        else:
            print("\n두 방법의 결과가 동일하지 않습니다.")
    except Exception as e:
        print(f"예상치 못한 오류가 발생했습니다: {e}")

        

if __name__ == "__main__":
    main()

