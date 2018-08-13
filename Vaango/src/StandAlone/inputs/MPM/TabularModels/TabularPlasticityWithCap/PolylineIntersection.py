def line_intersection(line1, line2):

    x11 = line1[0][0]
    y11 = line1[0][1]
    x12 = line1[1][0]
    y12 = line1[1][1]

    x21 = line2[0][0]
    y21 = line2[0][1]
    x22 = line2[1][0]
    y22 = line2[1][1]

    a11 = x12 - x11
    a12 = -(x22 - x21)
    a21 = y12 - y11
    a22 = -(y22 - y21)

    b1 = x21 - x11
    b2 = y21 - y11

    detA = a11*a22 - a12*a21
    if detA == 0:
       return None

    t1 = (a22*b1 - a12*b2)/detA
    t2 = (-a21*b1 + a11*b2)/detA

    x1 = (1-t1)*x11 + t1*x12
    y1 = (1-t1)*y11 + t1*y12

    #x2 = (1-t2)*x21 + t2*x22
    #y2 = (1-t2)*y21 + t2*y22

    #if (t1 >= 0 and t1 <= 1):
    #  print('t1 = ', t1, 'x1 = ', x1, 'y1 = ', y1)
    #  print('t2 = ', t2, 'x2 = ', x2, 'y2 = ', y2)
    #if (t2 >= 0 and t2 <= 1):
    #  print('t1 = ', t1, 'x1 = ', x1, 'y1 = ', y1)
    #  print('t2 = ', t2, 'x2 = ', x2, 'y2 = ', y2)

    return (x1, y1, t1, t2)

def poly_intersection(poly1, poly2):

    for i, p1_first_point in enumerate(poly1[:-1]):
        p1_second_point = poly1[i + 1]

        for j, p2_first_point in enumerate(poly2[:-1]):
            p2_second_point = poly2[j + 1]

            pt = line_intersection((p1_first_point, p1_second_point), (p2_first_point, p2_second_point))
            t1 = (pt[0] - p2_first_point[0])/(p2_second_point[0] - p2_first_point[0])
            t2 = (pt[1] - p2_first_point[1])/(p2_second_point[1] - p2_first_point[1])
            if (pt):
              if (pt[2] >= 0 and pt[2] <= 1 and pt[3] >= 0 and pt[3] <= 1):
                print(pt)
                #print(j, p1_first_point, p1_second_point, p2_first_point, p2_second_point, pt, t)
                return [True, pt, j]

    return False
