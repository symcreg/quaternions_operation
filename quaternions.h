//
// Created by symc on 2024/6/18.
//

#define FAST_INV_SQRT

#ifndef QUATERNIONS_QUATERNIONS_H
#define QUATERNIONS_QUATERNIONS_H
#include <math.h>
typedef struct {
    float x;
    float y;
    float z;
    float w;
} Quaternion;
void quaternion_add(Quaternion *result, Quaternion *a, Quaternion *b);
void quaternion_sub(Quaternion *result, Quaternion *a, Quaternion *b);
void quaternion_mul(Quaternion *result, Quaternion *a, Quaternion *b);
void quaternion_div(Quaternion *result, Quaternion *a, Quaternion *b);
void quaternion_conjugate(Quaternion *result, Quaternion *a);
void quaternion_inverse(Quaternion *result, Quaternion *a);
void quaternion_normalize(Quaternion *result, Quaternion *a);
void quaternion_rotate(Quaternion *result, Quaternion *a, Quaternion *b);
void quaternion_from_axis_angle(Quaternion *result, float x, float y, float z, float angle);
void quaternion_to_axis_angle(Quaternion *a, float *x, float *y, float *z, float *angle);
void quaternion_from_euler(Quaternion *result, float x, float y, float z);
void quaternion_to_euler(Quaternion *a, float *x, float *y, float *z);
void quaternion_from_matrix(Quaternion *result, float *matrix);
void quaternion_to_matrix(Quaternion *a, float *matrix);
void quaternion_slerp(Quaternion *result, Quaternion *a, Quaternion *b, float t);
void quaternion_lerp(Quaternion *result, Quaternion *a, Quaternion *b, float t);
void quaternion_nlerp(Quaternion *result, Quaternion *a, Quaternion *b, float t);
void quaternion_look_rotation(Quaternion *result, float x, float y, float z, float up_x, float up_y, float up_z);
void quaternion_update_from_gyroscope(Quaternion *result, Quaternion *q, float gx, float gy, float gz, float delta_t);
static inline float fastInvSqrt(float x) {
    float x2 = x * 0.5f;
    int i = *(int *) &x;
    i = 0x5f3759df - (i >> 1);
    x = *(float *) &i;
    x = x * (1.5f - (x2 * x * x));
    return x;
}
static inline float fastSqrt(float x) {
    return x * fastInvSqrt(x);
}
void quaternion_add(Quaternion *result, Quaternion *a, Quaternion *b) {
    result->x = a->x + b->x;
    result->y = a->y + b->y;
    result->z = a->z + b->z;
    result->w = a->w + b->w;
}
void quaternion_sub(Quaternion *result, Quaternion *a, Quaternion *b) {
    result->x = a->x - b->x;
    result->y = a->y - b->y;
    result->z = a->z - b->z;
    result->w = a->w - b->w;
}
void quaternion_mul(Quaternion *result, Quaternion *a, Quaternion *b) {
    result->x = a->w * b->x + a->x * b->w + a->y * b->z - a->z * b->y;
    result->y = a->w * b->y - a->x * b->z + a->y * b->w + a->z * b->x;
    result->z = a->w * b->z + a->x * b->y - a->y * b->x + a->z * b->w;
    result->w = a->w * b->w - a->x * b->x - a->y * b->y - a->z * b->z;
}
void quaternion_div(Quaternion *result, Quaternion *a, Quaternion *b) {
    Quaternion inv;
    quaternion_inverse(&inv, b);
    quaternion_mul(result, a, &inv);
}
void quaternion_conjugate(Quaternion *result, Quaternion *a) {
    result->x = -a->x;
    result->y = -a->y;
    result->z = -a->z;
    result->w = a->w;
}
void quaternion_inverse(Quaternion *result, Quaternion *a) {
    Quaternion conj;
    quaternion_conjugate(&conj, a);
    quaternion_normalize(result, &conj);
}
void quaternion_normalize(Quaternion *result, Quaternion *a) {
    float norm = a->x * a->x + a->y * a->y + a->z * a->z + a->w * a->w;
#ifdef FAST_INV_SQRT
    float inv_norm = fastInvSqrt(norm);
    result->x *= inv_norm;
    result->y *= inv_norm;
    result->z *= inv_norm;
    result->w *= inv_norm;
#else
    float inv_norm = 1.0f / sqrtf(norm);
    result->x *= inv_norm;
    result->y *= inv_norm;
    result->z *= inv_norm;
    result->w *= inv_norm;
#endif
}
void quaternion_rotate(Quaternion *result, Quaternion *a, Quaternion *b) {
    Quaternion q;
    Quaternion q_inv;
    quaternion_mul(&q, a, b);
    quaternion_inverse(&q_inv, a);
    quaternion_mul(result, &q, &q_inv);
}
void quaternion_from_axis_angle(Quaternion *result, float x, float y, float z, float angle) {
    float half_angle = angle * 0.5f;
    float s = sinf(half_angle);
    result->x = x * s;
    result->y = y * s;
    result->z = z * s;
    result->w = cosf(half_angle);
}
void quaternion_to_axis_angle(Quaternion *a, float *x, float *y, float *z, float *angle) {
#ifdef FAST_INV_SQRT
    float scale = fastSqrt(a->x * a->x + a->y * a->y + a->z * a->z);
#else
    float scale = sqrtf(a->x * a->x + a->y * a->y + a->z * a->z);
#endif
    float inv_scale = 1.0f / scale;
    if (scale == 0.0f) {
        *x = 0.0f;
        *y = 0.0f;
        *z = 0.0f;
        *angle = 0.0f;
    } else {
        *x = a->x * inv_scale;
        *y = a->y * inv_scale;
        *z = a->z * inv_scale;
        *angle = 2.0f * acosf(a->w);
    }
}
void quaternion_from_euler(Quaternion *result, float x, float y, float z) {
    float half_x = x * 0.5f;
    float half_y = y * 0.5f;
    float half_z = z * 0.5f;
    float c1 = cosf(half_x);
    float c2 = cosf(half_y);
    float c3 = cosf(half_z);
    float s1 = sinf(half_x);
    float s2 = sinf(half_y);
    float s3 = sinf(half_z);
    result->x = s1 * c2 * c3 + c1 * s2 * s3;
    result->y = c1 * s2 * c3 - s1 * c2 * s3;
    result->z = c1 * c2 * s3 + s1 * s2 * c3;
    result->w = c1 * c2 * c3 - s1 * s2 * s3;
}
void quaternion_to_euler(Quaternion *a, float *x, float *y, float *z) {
    float sqw = a->w * a->w;
    float sqx = a->x * a->x;
    float sqy = a->y * a->y;
    float sqz = a->z * a->z;
    *x = atan2f(2.0f * (a->x * a->w - a->y * a->z), sqx - sqy - sqz + sqw);
    *y = asinf(2.0f * (a->x * a->z + a->y * a->w));
    *z = atan2f(2.0f * (a->z * a->w - a->x * a->y), -sqx - sqy + sqz + sqw);
}
void quaternion_from_matrix(Quaternion *result, float *matrix) {
    float trace = matrix[0] + matrix[5] + matrix[10];
    if (trace > 0.0f) {
        float s = 0.5f / sqrtf(trace + 1.0f);
        result->w = 0.25f / s;
        result->x = (matrix[9] - matrix[6]) * s;
        result->y = (matrix[2] - matrix[8]) * s;
        result->z = (matrix[4] - matrix[1]) * s;
    } else {
        if (matrix[0] > matrix[5] && matrix[0] > matrix[10]) {
            float s = 2.0f * sqrtf(1.0f + matrix[0] - matrix[5] - matrix[10]);
            result->w = (matrix[9] - matrix[6]) / s;
            result->x = 0.25f * s;
            result->y = (matrix[1] + matrix[4]) / s;
            result->z = (matrix[2] + matrix[8]) / s;
        } else if (matrix[5] > matrix[10]) {
            float s = 2.0f * sqrtf(1.0f + matrix[5] - matrix[0] - matrix[10]);
            result->w = (matrix[2] - matrix[8]) / s;
            result->x = (matrix[1] + matrix[4]) / s;
            result->y = 0.25f * s;
            result->z = (matrix[6] + matrix[9]) / s;
        } else {
            float s = 2.0f * sqrtf(1.0f + matrix[10] - matrix[0] - matrix[5]);
            result->w = (matrix[4] - matrix[1]) / s;
            result->x = (matrix[2] + matrix[8]) / s;
            result->y = (matrix[6] + matrix[9]) / s;
            result->z = 0.25f * s;
        }
    }
}
void quaternion_to_matrix(Quaternion *a, float *matrix) {
    float xx = a->x * a->x;
    float xy = a->x * a->y;
    float xz = a->x * a->z;
    float xw = a->x * a->w;
    float yy = a->y * a->y;
    float yz = a->y * a->z;
    float yw = a->y * a->w;
    float zz = a->z * a->z;
    float zw = a->z * a->w;
    matrix[0] = 1.0f - 2.0f * (yy + zz);
    matrix[1] = 2.0f * (xy - zw);
    matrix[2] = 2.0f * (xz + yw);
    matrix[3] = 0.0f;
    matrix[4] = 2.0f * (xy + zw);
    matrix[5] = 1.0f - 2.0f * (xx + zz);
    matrix[6] = 2.0f * (yz - xw);
    matrix[7] = 0.0f;
    matrix[8] = 2.0f * (xz - yw);
    matrix[9] = 2.0f * (yz + xw);
    matrix[10] = 1.0f - 2.0f * (xx + yy);
    matrix[11] = 0.0f;
    matrix[12] = 0.0f;
    matrix[13] = 0.0f;
    matrix[14] = 0.0f;
    matrix[15] = 1.0f;
}
void quaternion_slerp(Quaternion *result, Quaternion *a, Quaternion *b, float t) {
    float cos_half_theta = a->x * b->x + a->y * b->y + a->z * b->z + a->w * b->w;
    if (cos_half_theta < 0.0f) {
        a->x = -a->x;
        a->y = -a->y;
        a->z = -a->z;
        a->w = -a->w;
        cos_half_theta = -cos_half_theta;
    }
    if (cos_half_theta >= 1.0f) {
        result->x = a->x;
        result->y = a->y;
        result->z = a->z;
        result->w = a->w;
    } else {
        float half_theta = acosf(cos_half_theta);
        float sin_half_theta = sqrtf(1.0f - cos_half_theta * cos_half_theta);
        if (fabsf(sin_half_theta) < 0.001f) {
            result->x = a->x * 0.5f + b->x * 0.5f;
            result->y = a->y * 0.5f + b->y * 0.5f;
            result->z = a->z * 0.5f + b->z * 0.5f;
            result->w = a->w * 0.5f + b->w * 0.5f;
        } else {
            float ratio_a = sinf((1.0f - t) * half_theta) / sin_half_theta;
            float ratio_b = sinf(t * half_theta) / sin_half_theta;
            result->x = a->x * ratio_a + b->x * ratio_b;
            result->y = a->y * ratio_a + b->y * ratio_b;
            result->z = a->z * ratio_a + b->z * ratio_b;
            result->w = a->w * ratio_a + b->w * ratio_b;
        }
    }
}
void quaternion_lerp(Quaternion *result, Quaternion *a, Quaternion *b, float t) {
    result->x = a->x + (b->x - a->x) * t;
    result->y = a->y + (b->y - a->y) * t;
    result->z = a->z + (b->z - a->z) * t;
    result->w = a->w + (b->w - a->w) * t;
}
void quaternion_nlerp(Quaternion *result, Quaternion *a, Quaternion *b, float t) {
    Quaternion temp;
    quaternion_lerp(&temp, a, b, t);
    quaternion_normalize(result, &temp);
}
void quaternion_look_rotation(Quaternion *result, float x, float y, float z, float up_x, float up_y, float up_z) {
    float forward_x = x;
    float forward_y = y;
    float forward_z = z;
    float up_x_ = up_x;
    float up_y_ = up_y;
    float up_z_ = up_z;
#ifdef FAST_INV_SQRT
    float forward_length = fastSqrt(forward_x * forward_x + forward_y * forward_y + forward_z * forward_z);
    float up_length = fastSqrt(up_x_ * up_x_ + up_y_ * up_y_ + up_z_ * up_z_);
#else
    float forward_length = sqrtf(forward_x * forward_x + forward_y * forward_y + forward_z * forward_z);
    float up_length = sqrtf(up_x_ * up_x_ + up_y_ * up_y_ + up_z_ * up_z_);
#endif
    float inv_forward_length = 1.0f / forward_length;
    float inv_up_length = 1.0f / up_length;
    if (forward_length < 0.0001f) {
        forward_x = 0.0f;
        forward_y = 0.0f;
        forward_z = 1.0f;
    } else {
        forward_x *= inv_forward_length;
        forward_y *= inv_forward_length;
        forward_z *= inv_forward_length;
    }
    if (up_length < 0.0001f) {
        up_x_ = 0.0f;
        up_y_ = 1.0f;
        up_z_ = 0.0f;
    } else {
        up_x_ *= inv_up_length;
        up_y_ *= inv_up_length;
        up_z_ *= inv_up_length;
    }
    float right_x = up_y_ * forward_z - up_z_ * forward_y;
    float right_y = up_z_ * forward_x - up_x_ * forward_z;
    float right_z = up_x_ * forward_y - up_y_ * forward_x;
#ifdef FAST_INV_SQRT
    float right_length = fastSqrt(right_x * right_x + right_y * right_y + right_z * right_z);
#else
    float right_length = sqrtf(right_x * right_x + right_y * right_y + right_z * right_z);
#endif
    float inv_right_length = 1.0f / right_length;
    if (right_length < 0.0001f) {
        right_x = 1.0f;
        right_y = 0.0f;
        right_z = 0.0f;
    } else {
        right_x *= inv_right_length;
        right_y *= inv_right_length;
        right_z *= inv_right_length;
    }
    up_x_ = forward_y * right_z - forward_z * right_y;
    up_y_ = forward_z * right_x - forward_x * right_z;
    up_z_ = forward_x * right_y - forward_y * right_x;
    float m00 = right_x;
    float m01 = right_y;
    float m02 = right_z;
    float m10 = up_x_;
    float m11 = up_y_;
    float m12 = up_z_;
    float m20 = forward_x;
    float m21 = forward_y;
    float m22 = forward_z;
    float num8 = m00 + m11 + m22;
    if (num8 > 0.0f) {
#ifdef FAST_INV_SQRT
        float num = fastSqrt(num8 + 1.0f);
#else
        float num = sqrtf(num8 + 1.0f);
#endif
        result->w = num * 0.5f;
        num = 0.5f / num;
        result->x = (m12 - m21) * num;
        result->y = (m20 - m02) * num;
        result->z = (m01 - m10) * num;
    } else if (m00 >= m11 && m00 >= m22) {
#ifdef FAST_INV_SQRT
        float num7 = fastSqrt(1.0f + m00 - m11 - m22);
#else
        float num7 = sqrtf(1.0f + m00 - m11 - m22);
#endif
        float num4 = 0.5f / num7;
        result->x = 0.5f * num7;
        result->y = (m01 + m10) * num4;
        result->z = (m02 + m20) * num4;
        result->w = (m12 - m21) * num4;
    } else if (m11 > m22) {
#ifdef FAST_INV_SQRT
        float num6 = fastSqrt(1.0f + m11 - m00 - m22);
#else
        float num6 = sqrtf(1.0f + m11 - m00 - m22);
#endif
        float num3 = 0.5f / num6;
        result->x = (m10 + m01) * num3;
        result->y = 0.5f * num6;
        result->z = (m21 + m12) * num3;
        result->w = (m20 - m02) * num3;
    } else {
#ifdef FAST_INV_SQRT
        float num5 = fastSqrt(1.0f + m22 - m00 - m11);
#else
        float num5 = sqrtf(1.0f + m22 - m00 - m11);
#endif
        float num2 = 0.5f / num5;
        result->x = (m20 + m02) * num2;
        result->y = (m21 + m12) * num2;
        result->z = 0.5f * num5;
        result->w = (m01 - m10) * num2;
    }
}
void quaternion_update_from_gyroscope(Quaternion *result, Quaternion *q, float gx, float gy, float gz, float delta_t) {
    float half_t = delta_t * 0.5f;
    float sin_half_t = sinf(half_t);
    float cos_half_t = cosf(half_t);
    Quaternion delta_q;
    delta_q.x = gx * sin_half_t;
    delta_q.y = gy * sin_half_t;
    delta_q.z = gz * sin_half_t;
    delta_q.w = cos_half_t;
    quaternion_mul(result, q, &delta_q);
}
#endif //QUATERNIONS_QUATERNIONS_H
