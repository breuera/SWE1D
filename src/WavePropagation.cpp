/**
 * @file
 *  This file is part of SWE1D
 *
 *  SWE1D is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  SWE1D is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with SWE1D.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Diese Datei ist Teil von SWE1D.
 *
 *  SWE1D ist Freie Software: Sie koennen es unter den Bedingungen
 *  der GNU General Public License, wie von der Free Software Foundation,
 *  Version 3 der Lizenz oder (nach Ihrer Option) jeder spaeteren
 *  veroeffentlichten Version, weiterverbreiten und/oder modifizieren.
 *
 *  SWE1D wird in der Hoffnung, dass es nuetzlich sein wird, aber
 *  OHNE JEDE GEWAEHELEISTUNG, bereitgestellt; sogar ohne die implizite
 *  Gewaehrleistung der MARKTFAEHIGKEIT oder EIGNUNG FUER EINEN BESTIMMTEN
 *  ZWECK. Siehe die GNU General Public License fuer weitere Details.
 *
 *  Sie sollten eine Kopie der GNU General Public License zusammen mit diesem
 *  Programm erhalten haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>.
 * 
 * @copyright 2013 Technische Universitaet Muenchen
 * @author Sebastian Rettenberger <rettenbs@in.tum.de>
 */

#include "WavePropagation.h"

T WavePropagation::computeNumericalFluxes()
{
	float maxWaveSpeed = 0.f;

	// Loop over all edges
	for (unsigned int i = 1; i < m_size+2; i++) {
		T maxEdgeSpeed;

    // convert from soa t aos
    double l_variablesLeft[3];
    double l_variablesRight[3];
    double l_netUpdatesLeft[3];
    double l_netUpdatesRight[3];
    double l_waveSpeeds[3];

    l_variablesLeft[0] = m_h[i-1];
    l_variablesLeft[1] = m_hu[i-1];
    l_variablesLeft[2] = 0;

    l_variablesRight[0] = m_h[i];
    l_variablesRight[1] = m_hu[i];
    l_variablesRight[2] = 0;

    // call the GeoClaw solver
    c_bind_geoclaw_riemann_aug_JCP( 1,
                                    l_variablesLeft, l_variablesRight,
                                    0.0001, 9.81,
                                    l_netUpdatesLeft,l_netUpdatesRight,
                                    l_waveSpeeds );

    // convert from aos to soa
    m_hNetUpdatesLeft[i-1]  = l_netUpdatesLeft[0];
    m_huNetUpdatesLeft[i-1] = l_netUpdatesLeft[1];

    m_hNetUpdatesRight[i-1]  = l_netUpdatesRight[0];
    m_huNetUpdatesRight[i-1] = l_netUpdatesRight[1];

    // compute maximum wave speed from first and third wave family
    maxEdgeSpeed = std::max( std::abs( l_waveSpeeds[0] ), std::abs( l_waveSpeeds[2] ) );

		// Update maxWaveSpeed
		if (maxEdgeSpeed > maxWaveSpeed)
			maxWaveSpeed = maxEdgeSpeed;
	}

	// Compute CFL condition
	T maxTimeStep = m_cellSize/maxWaveSpeed * .4f;

	return maxTimeStep;
}

void WavePropagation::updateUnknowns(T dt)
{
	// Loop over all inner cells
	for (unsigned int i = 1; i < m_size+1; i++) {
        m_h[i] -=  dt/m_cellSize * (m_hNetUpdatesRight[i-1] + m_hNetUpdatesLeft[i]);
        m_hu[i] -= dt/m_cellSize * (m_huNetUpdatesRight[i-1] + m_huNetUpdatesLeft[i]);
	}
}

void WavePropagation::setOutflowBoundaryConditions()
{
	m_h[0] = m_h[1]; m_h[m_size+1] = m_h[m_size];
	m_hu[0] = m_hu[1]; m_hu[m_size+1] = m_hu[m_size];
}
