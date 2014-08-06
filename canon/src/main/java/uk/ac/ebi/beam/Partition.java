/*
 * Copyright (c) 2014, European Bioinformatics Institute (EMBL-EBI)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met: 
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer. 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution. 
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies, 
 * either expressed or implied, of the FreeBSD Project.
 */

package uk.ac.ebi.beam;

import java.util.Arrays;

/**
 * @author John May
 */
final class Partition {

    final int[] size;
    final int[] order;

    Partition(int[] order, int[] size) {
        this.order = order;
        this.size  = size;
    }
    
    boolean refine(Comparator comparator) {
        // split cells using comparator
        
        return false;    
    }

    @Override public String toString() {
        StringBuilder sb = new StringBuilder();
        int i = 0;
        sb.append('[');
        while (i < order.length) {
            if (i > 0) sb.append('|');
            int n = size[i];
            while (n-- > 0) {
                sb.append(order[i++]);
                if (n > 0)
                    sb.append(' ');
            }    
        }
        sb.append(']');
        return sb.toString();
    }

    static Partition copyOf(Partition org) {
        int n = org.order.length;
        return new Partition(Arrays.copyOf(org.order, n),
                             Arrays.copyOf(org.size, n));
    }

    static Partition unit(int n) {
        int[] order = new int[n];
        int[] size  = new int[n];
        
        for (int i = 0; i < n; i++)
            order[i] = i;
        
        size[0] = n;
        
        return new Partition(order, size);
    }
}
