#include "me_estimator.h"
#include "metric.h"
#include <iostream>
#include <array>

MotionEstimator::MotionEstimator(size_t width, size_t height, unsigned char quality, bool use_half_pixel)
: width(width)
, height(height)
, quality(quality)
, use_half_pixel(use_half_pixel)
, width_ext(width + 2 * BORDER)
, num_blocks_hor((width + BLOCK_SIZE - 1) / BLOCK_SIZE)
, num_blocks_vert((height + BLOCK_SIZE - 1) / BLOCK_SIZE)
, first_row_offset(width_ext * BORDER + BORDER)
, me_field(num_blocks_hor, num_blocks_vert, BLOCK_SIZE)
, width_borders(width + 2 * BORDER)
, height_borders(height + 2 * BORDER)
{
    cur_Y_borders = new unsigned char[width_borders * height_borders];
    prev_Y_borders = new unsigned char[width_borders * height_borders];
    if (use_half_pixel) {
        prev_Y_up_borders = new unsigned char[width_borders * height_borders];
        prev_Y_up_left_borders = new unsigned char[width_borders * height_borders];
        prev_Y_left_borders = new unsigned char[width_borders * height_borders];
    }
    
#if 0    
    std::cout << this->width << " " << this->height  << " " << this->quality << "\n";
    std::cout << this->use_half_pixel << " " << this->width_ext << " ";
    std::cout << this->num_blocks_hor << "\n";
    std::cout << this->num_blocks_vert << " " << this->first_row_offset << "\n";
    std::cout << this->width_borders << " " << height_borders << "\n";
#endif
    
}

MotionEstimator::~MotionEstimator() {
    delete[] cur_Y_borders;
    delete[] prev_Y_borders;
    if (use_half_pixel) {
        delete[] prev_Y_up_borders;
        delete[] prev_Y_up_left_borders;
        delete[] prev_Y_left_borders;
    }
}

void MotionEstimator::CEstimate(const unsigned char* cur_Y,
                                const unsigned char* prev_Y,
                                const uint8_t* prev_Y_up,
                                const uint8_t* prev_Y_left,
                                const uint8_t* prev_Y_upleft,
                                MEField& mvectors) {
    
    std::unordered_map<ShiftDir, const uint8_t*, ShiftDirHash> prev_map {
        { ShiftDir::NONE, prev_Y }
    };
    
    if (use_half_pixel) {
        prev_map.emplace(ShiftDir::UP, prev_Y_up);
        prev_map.emplace(ShiftDir::LEFT, prev_Y_left);
        prev_map.emplace(ShiftDir::UPLEFT, prev_Y_upleft);
    }
    
    const int    sad_thrashold  = (this->quality == 100) ? 100
        : (this->quality == 80) ? 100
        : (this->quality == 60) ? 200
        : (this->quality == 40) ? 300 : 500;
    const double disp_thrashold = (this->quality == 100) ? 350
        : (this->quality == 80) ? 400
        : (this->quality == 60) ? 425
        : (this->quality == 40) ? 450 : 475;
    
    
    static int frameID = -1;
    ++frameID;
    
    for (size_t i = 0; i < num_blocks_vert; ++i) {
        for (size_t j = 0; j < num_blocks_hor; ++j) {
            
            const auto hor_offset = j * BLOCK_SIZE;
            const auto vert_offset = first_row_offset + i * BLOCK_SIZE * width_ext;
            auto cur = cur_Y + vert_offset + hor_offset;
            auto prev = prev_Y + vert_offset + hor_offset;
            
            MV best_vector;
            
            // best_vector x, y are 0
            best_vector.error = GetErrorSAD(cur, prev, width_ext, 16);
            
            if (best_vector.error > sad_thrashold) {
                
                if (GetDispersion(cur, width_ext, 16) > disp_thrashold) {
                    const int &l_width = width_ext;
                    
                    // lambda
                    auto ds_check_new = [&prev, &cur, &l_width](MV        &l_vector,
                                                                const int  l_block_size,
                                                                const int  l_len_coords,
                                                                const int *l_xs,
                                                                const int *l_ys
                                                                ) -> int {
                        
                        int best_move = -1;
                        
                        for (int ii = 0; ii < l_len_coords; ++ii) {
                            const auto comp = prev
                                + (l_vector.y + l_ys[ii]) * l_width
                                + (l_vector.x + l_xs[ii]);
                            
                            const int error = GetErrorSAD(cur, comp, l_width, l_block_size);
                            
                            if (error < l_vector.error) {
                                best_move = ii;
                                l_vector.error = error;
                            }
                        }
                        
                        return best_move;
                    };
                    
                    // lambda
                    auto ds_full8 = [&prev, &cur, &l_width] (MV            &l_vector,
                                                             const int      l_block_size,
                                                             const int      l_ds_step,
                                                             const uint8_t *center,
                                                             const int      l_center_x = 0,
                                                             const int      l_center_y = 0
                                                             ) -> void {
                        const int xs[] = {
                            -2 * l_ds_step,
                            -l_ds_step, -l_ds_step,
                            0, 0,
                            l_ds_step, l_ds_step,
                            2 * l_ds_step
                        };
                        
                        const int ys[] = {
                            0,
                            -l_ds_step, l_ds_step,
                            -2 * l_ds_step, 2 * l_ds_step,
                            -l_ds_step, l_ds_step,
                            0
                        };
                        
                        for (int ii = 0; ii < 8; ++ii) {
                            const auto x = xs[ii];
                            const auto y = ys[ii];
                            const auto comp = center + y * l_width + x;
                            
                            const int error = GetErrorSAD(cur, comp, l_width, l_block_size);
                            
                            if (error < l_vector.error) {
                                l_vector.x = l_center_x + x;
                                l_vector.y = l_center_y + y;
                                l_vector.error = error;
                            }
                        }
                    };
                    
                    const int lds_step = 2;
                    ds_full8(best_vector, BLOCK_SIZE, lds_step, prev);
                    
                    auto b_x = best_vector.x;
                    auto b_y = best_vector.y;
                    
                    while (b_x || b_y) {
                        if ( abs(b_x) == lds_step ) {
                            const int xs[] = {
                                b_x, 2 * b_x, 0
                            };
                            const int ys[] = {
                                b_y, 0, 2 * b_y
                            };
                            const int best_move = ds_check_new(best_vector, BLOCK_SIZE, 3, xs, ys);
                            
                            if (best_move < 0) {
                                break;
                            } else {
                                b_x = xs[best_move];
                                b_y = ys[best_move];
                            }
                            
                        } else {
                            const int xs[] = {
                                b_x, b_y, -b_y, (b_x + b_y) / 2,
                                ((!b_x) ? -1 : 1) * (b_x + b_y) / 2
                            };
                            const int ys[] = {
                                b_y, b_x, -b_x, (b_x + b_y) / 2,
                                ((!b_y) ? -1 : 1) * (b_y + b_x) / 2
                            };
                            const int best_move = ds_check_new(best_vector, BLOCK_SIZE, 5, xs, ys);
                            
                            if (best_move < 0) {
                                break;
                            } else {
                                b_x = xs[best_move];
                                b_y = ys[best_move];
                            }
                        }
                        
                        best_vector.x += b_x;
                        best_vector.y += b_y;
                    }
                    
                    // Refine vector with SDS
                    const auto rough_vector = prev + best_vector.y * width_ext + best_vector.x;
                    ds_full8(best_vector, BLOCK_SIZE, lds_step / 2, rough_vector, best_vector.x, best_vector.y);
                    
                    if (frameID == 20) {
                        std::cout << "found vector for frame #" << frameID << "; ";
                        std::cout << "(" << i << "; " << j << ") " <<  best_vector.x << " " << best_vector.y << "\n";
                    }
                    
                    // Split into four subvectors if the error is too large
                    if (1) {
                        best_vector.Split();
                        
                        for (int h = 0; h < 4; ++h) {
                            auto& subvector = best_vector.SubVector(h);
                            
                            const auto hor_offset  = j * BLOCK_SIZE + ((h & 1) ? BLOCK_SIZE / 2 : 0);
                            const auto vert_offset = first_row_offset
                                + (i * BLOCK_SIZE + ((h > 1) ? BLOCK_SIZE / 2 : 0)) * width_ext;
                            cur  = cur_Y  + vert_offset + hor_offset;
                            prev = prev_Y + vert_offset + hor_offset;
                            
                            subvector.error = GetErrorSAD(cur, prev, width_ext, 8);
                            
                            if (subvector.error > 0) {
                                if (GetDispersion(cur, width_ext, 8) > 0) {
                                    ds_full8(subvector, 8, lds_step, prev);
                                    auto b_x = subvector.x;
                                    auto b_y = subvector.y;
                                    
                                    while (b_x || b_y) {
                                        if ( abs(b_x) == lds_step ) {
                                            const int xs[] = {
                                                b_x, 2 * b_x, 0
                                            };
                                            const int ys[] = {
                                                b_y, 0, 2 * b_y
                                            };
                                            const int best_move = ds_check_new(subvector, 8, 3, xs, ys);
                                            
                                            if (best_move < 0) {
                                                break;
                                            } else {
                                                b_x = xs[best_move];
                                                b_y = ys[best_move];
                                            }
                                            
                                        } else {
                                            const int xs[] = {
                                                b_x, b_y, -b_y, (b_x + b_y) / 2,
                                                ((!b_x) ? -1 : 1) * (b_x + b_y) / 2
                                            };
                                            const int ys[] = {
                                                b_y, b_x, -b_x, (b_x + b_y) / 2,
                                                ((!b_y) ? -1 : 1) * (b_y + b_x) / 2
                                            };
                                            const int best_move = ds_check_new(subvector, 8, 5, xs, ys);
                                            
                                            if (best_move < 0) {
                                                break;
                                            } else {
                                                b_x = xs[best_move];
                                                b_y = ys[best_move];
                                            }
                                        }
                                        
                                        subvector.x += b_x;
                                        subvector.y += b_y;
                                    }
                                    
                                    
                                    // Refine vector with SDS
                                    const auto rough_vector = prev + subvector.y * width_ext + subvector.x;
                                    ds_full8(subvector,
                                             BLOCK_SIZE,
                                             lds_step / 2,
                                             rough_vector,
                                             subvector.x,
                                             subvector.y);
                                }
                            }
                        }
                        
                        if (best_vector.SubVector(0).error
                            + best_vector.SubVector(1).error
                            + best_vector.SubVector(2).error
                            + best_vector.SubVector(3).error >= best_vector.error
                            ) {
                            best_vector.Unsplit();
                        }
                    }
                    
                    mvectors.set_mv(j, i, best_vector);
                }
            }
        }
    }
}

void MotionEstimator::BruteForce(const unsigned char* cur_Y,
                                 const unsigned char* prev_Y,
                                 const uint8_t* prev_Y_up,
                                 const uint8_t* prev_Y_left,
                                 const uint8_t* prev_Y_upleft,
                                 MEField& mvectors) {
    std::unordered_map<ShiftDir, const uint8_t*, ShiftDirHash> prev_map {
        { ShiftDir::NONE, prev_Y }
    };
    
    if (use_half_pixel) {
        prev_map.emplace(ShiftDir::UP, prev_Y_up);
        prev_map.emplace(ShiftDir::LEFT, prev_Y_left);
        prev_map.emplace(ShiftDir::UPLEFT, prev_Y_upleft);
    }
    
    for (size_t i = 0; i < num_blocks_vert; ++i) {
        for (size_t j = 0; j < num_blocks_hor; ++j) {
            const auto hor_offset = j * BLOCK_SIZE;
            const auto vert_offset = first_row_offset + i * BLOCK_SIZE * width_ext;
            const auto cur = cur_Y + vert_offset + hor_offset;
            
            MV best_vector;
            best_vector.error = std::numeric_limits<long>::max();
            
            // Brute force
            for (const auto& prev_pair : prev_map) {
                const auto prev = prev_pair.second + vert_offset + hor_offset;
                for (int y = -BORDER; y <= BORDER; ++y) {
                    for (int x = -BORDER; x <= BORDER; ++x) {
                        const auto comp = prev + y * width_ext + x;
                        const int error = GetErrorSAD_16x16(cur, comp, width_ext);
                        if (error < best_vector.error) {
                            best_vector.x = x;
                            best_vector.y = y;
                            best_vector.shift_dir = prev_pair.first;
                            best_vector.error = error;
                        }
                    }
                }
            }
            
            // Split into four subvectors if the error is too large
            if (best_vector.error > 1000) {
                best_vector.Split();
                
                for (int h = 0; h < 4; ++h) {
                    auto& subvector = best_vector.SubVector(h);
                    subvector.error = std::numeric_limits<long>::max();
                    
                    const auto hor_offset = j * BLOCK_SIZE + ((h & 1) ? BLOCK_SIZE / 2 : 0);
                    const auto vert_offset = first_row_offset +
                        (i * BLOCK_SIZE + ((h > 1) ? BLOCK_SIZE / 2 : 0)) * width_ext;
                    const auto cur = cur_Y + vert_offset + hor_offset;
                    
                    for (const auto& prev_pair : prev_map) {
                        const auto prev = prev_pair.second + vert_offset + hor_offset;
                        
                        for (int y = -BORDER; y <= BORDER; ++y) {
                            for (int x = -BORDER; x <= BORDER; ++x) {
                                const auto comp = prev + y * width_ext + x;
                                const int error = GetErrorSAD_8x8(cur, comp, width_ext);
                                
                                if (error < subvector.error) {
                                    subvector.x = x;
                                    subvector.y = y;
                                    subvector.shift_dir = prev_pair.first;
                                    subvector.error = error;
                                }
                            }
                        }
                    }
                }
                
                if (best_vector.SubVector(0).error
                    + best_vector.SubVector(1).error
                    + best_vector.SubVector(2).error
                    + best_vector.SubVector(3).error > best_vector.error * 0.7
                    ) {
                    best_vector.Unsplit();
                }
            }
            mvectors.set_mv(j, i, best_vector);
        }
    }
}

void extend_with_borders(
                         unsigned char *input,
                         unsigned char *output,
                         size_t height,
                         size_t width,
                         size_t border_size
                         ) {
    // Copy frame to center of new 
    size_t new_width = width + 2 * border_size;
    auto p_output = output + new_width * border_size + border_size;
    auto p_input = input;
    for (size_t y = 0; y < height; ++y, p_output += 2 * border_size) {
        for (size_t x = 0; x < width; ++x, ++p_output, ++p_input) {
            *p_output = *p_input;
        }
    }
    
    // Left and right borders.
    p_output = output + new_width * border_size;
    for (size_t y = 0; y < height; ++y) {
        memset(p_output, p_output[border_size], border_size);
        p_output += border_size + width;
        memset(p_output, p_output[-1], border_size);
        p_output += border_size;
    }
    
    // Top and bottom borders.
    p_output = output;
    auto p_output_reference_row = p_output + new_width * border_size;
    
    for (size_t y = 0; y < border_size; ++y) {
        memcpy(p_output, p_output_reference_row, new_width);
        p_output += new_width;
    }
    p_output = output + new_width * (height + border_size);
    p_output_reference_row = p_output_reference_row - new_width;
    
    for (size_t y = 0; y < border_size; ++y) {
        memcpy(p_output, p_output_reference_row, new_width);
        p_output += new_width;
    }
}

void generate_subpixel_arrays(
                              unsigned char* input,
                              unsigned char* output_up,
                              unsigned char* output_left,
                              unsigned char* output_up_left,
                              size_t height,
                              size_t width
                              ) {
    for (size_t y = 0; y < height; ++y) {
        for (size_t x = 0; x < width; ++x) {
            size_t cur_pixel_pos = y * width + x;
            size_t left_pixel_pos = y * width + x - 1;
            size_t left_up_pixel_pos = (y - 1) * width + x - 1;
            size_t up_pixel_pos = (y - 1) * width + x;
            
            if (y > 0) {
                output_up[cur_pixel_pos] = (int(input[cur_pixel_pos]) + input[up_pixel_pos]) / 2;
            } else {
                output_up[cur_pixel_pos] = input[cur_pixel_pos];
            }
            if (x > 0) {
                output_left[cur_pixel_pos] = (int(input[cur_pixel_pos]) + input[left_pixel_pos]) / 2;
            } else {
                output_left[cur_pixel_pos] = input[cur_pixel_pos];
            }
            
            if (x > 0 && y > 0) {
                output_up_left[cur_pixel_pos] = (
                                                 int(input[cur_pixel_pos]) + 
                                                 input[left_pixel_pos] + 
                                                 input[left_up_pixel_pos] + 
                                                 input[up_pixel_pos]
                                                 ) / 4;
            } else if (y == 0) {
                output_up_left[cur_pixel_pos] = output_left[cur_pixel_pos];
            } else {
                output_up_left[cur_pixel_pos] = output_up[cur_pixel_pos];
            }
        }
    }   
}

MEField MotionEstimator::Estimate(
                                  py::array_t<unsigned char> cur_Y,
                                  py::array_t<unsigned char> prev_Y
                                  ) {
    
    extend_with_borders((unsigned char *)cur_Y.request().ptr,
                        cur_Y_borders, height, width, BORDER);
    extend_with_borders((unsigned char *)prev_Y.request().ptr,
                        prev_Y_borders, height, width, BORDER);
    
    if (cur_Y.size() != prev_Y.size()) {
        throw std::runtime_error("Input shapes must match");
    }
    
    if (use_half_pixel) {
        generate_subpixel_arrays(
                                 prev_Y_borders,
                                 prev_Y_up_borders,
                                 prev_Y_left_borders,
                                 prev_Y_up_left_borders,
                                 width_borders,
                                 height_borders
                                 );
    }
      
    MotionEstimator::BruteForce(
                                cur_Y_borders,
                                prev_Y_borders,
                                prev_Y_up_borders,
                                prev_Y_left_borders,
                                prev_Y_up_left_borders,
                                me_field
                                );
#if 0
    
    MotionEstimator::CEstimate(
                               cur_Y_borders,
                               prev_Y_borders,
                               prev_Y_up_borders,
                               prev_Y_left_borders,
                               prev_Y_up_left_borders,
                               me_field
                               );
    #endif

    return me_field;
    
}
