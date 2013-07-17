/*
 * Copyright (c) 2013, European Bioinformatics Institute (EMBL-EBI)
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

package uk.ac.ebi.grins;

/**
 * Defines the relative topology around a vertex (atom).
 *
 * @author John May
 */
abstract class Topology {

    /**
     * The vertex/atom which this topology describes.
     *
     * @return vertex
     * @throws IllegalArgumentException unknown topology
     */
    abstract int atom();

    /**
     * The configuration of the topology.
     *
     * @return configuration for this topology
     */
    abstract Configuration configuration();

    /**
     * Arrange the topology relative to a given ranking of vertices.
     *
     * @param rank ordering of vertices
     * @return a new topology with the neighbors arranged by the given rank
     */
    abstract Topology orderBy(int[] rank);

    /**
     * Transform the topology to one with the given {@literal mapping}.
     *
     * @param mapping the mapping used to transform the topology
     * @return a new topology with it's vertices mapped
     */
    abstract Topology transform(int[] mapping);

    /**
     * Specify unknown configuration on atom - there is no vertex data stored.
     *
     * @return unknown topology
     */
    static Topology unknown() {
        return UNKNOWN;
    }

    private static Topology UNKNOWN = new Topology() {
        @Override int atom() {
            throw new IllegalArgumentException("unknown topology");
        }

        @Override Configuration configuration() {
            return Configuration.UNKNOWN;
        }

        @Override Topology orderBy(int[] rank) {
            return this;
        }

        @Override Topology transform(int[] mapping) {
            return this;
        }
    };
}
