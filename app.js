// Constants
const P_ATMOSPHERE = 101325; // Pa
const D_ATMOSPHERE = 1.225; //standard atmospheric density
const V_BOTTLE     = 0.001;  // m³ (1 liter)
let   V_WATER      = 0.00015; // m³ (initial water volume)
const GAMMA        = 1.4;    // Adiabatic index
const RHO_WATER    = 1000;   // kg/m³
const NOZZLE_RADIUS= 0.0045; // m
const NOZZLE_AREA  = Math.PI * NOZZLE_RADIUS ** 2; // m²
const CD = 1;
const BOTTLE_MASS  = 0.15; //kg
const TEMP_INITIAL = 288; //in Kelvin
const GAS_CONSTANT = 287; // in Joules


function simulate(P_psi) {
    // Convert PSI to Pascals (gauge => absolute means + atmospheric,
    // but if your P_psi is "gauge" above atmosphere, just multiply by 6894.76)
    const P_INITIAL   = P_psi * 6894.76; 
    // Initial air volume
    const V_GAS_INITIAL = V_BOTTLE - V_WATER; 
    
    // We'll track volumes, pressures, times for plotting
    let volumes    = [];
    let pressures  = [];
    let times      = [];
    let areas      = [];
    let thrusts    = [];
    let airthrusts = [];
    let airtime = [];
    let temperatures = [];
    let totalArea  = 0;
    let time       = 0;
    let totalThrust = 0;

    // Current water volume & derived current air volume
    let V_water = V_WATER;          
    let V_gas   = V_GAS_INITIAL;    
    
    // We'll use a small step in water volume each iteration
    const dV_water = 1e-7;   // m³ (adjust as needed for speed/accuracy)
    
    // Keep track of "previous" pressure for trapezoidal integration
    let P_previous = P_INITIAL;
    let T_previous = TEMP_INITIAL;

    // Loop until the water is gone or the pressure hits atmosphere
    while (V_water > 0) {
        // Current air volume
        V_gas = V_BOTTLE - V_water;
        
        // Compute the current pressure using adiabatic expansion
        let P = P_INITIAL * Math.pow(V_GAS_INITIAL / V_gas, GAMMA);

        // If the pressure drops to atmospheric, rocket no longer expels water
        //if (P <= P_ATMOSPHERE) break;

        // Calculate water exit velocity using Bernoulli
        let velocity = CD * Math.sqrt((2 * (P)) / RHO_WATER);
        
        // Flow rate of water (m³/s)
        let flowRate = NOZZLE_AREA * velocity;
        let massFlowRate = flowRate * RHO_WATER;

        //Need to remove the force due to gravity from the thrust
        let waterMass = V_water * RHO_WATER;
        let gravityForce = (BOTTLE_MASS + waterMass) * 9.81;

        // calculate thrust
        let thrust = (massFlowRate * velocity) - gravityForce;

        // If the flow rate is effectively zero (just in case), break
        if (flowRate < 1e-12) break;

        // Time step for removing dV_water
        let dt = dV_water / flowRate;
        time += dt;
        
        //add up the thrust values
        totalThrust += thrust;

        // Update water volume (expelled some portion)
        V_water -= dV_water;
        if (V_water < 0) V_water = 0; // clamp

        //calculate the temperature at this point
        let new_gas = V_BOTTLE - V_water
        let temperature = T_previous * ((V_gas/(new_gas)) ** 0.4)

        // Compute area under P-V curve for this step 
        // (We approximate using the "air volume" change => dV_air = dV_water,
        //  but strictly the area is for the P vs AIR-volume curve.)
        // You could do something fancier to map dV_water -> dV_air, 
        // but for short steps this is close enough:
        let areaStep = ((P + P_previous) / 2) * (dV_water); 
        areas.push(areaStep);
        totalArea += areaStep;

        // Save data for plotting
        volumes.push(V_gas);   // we track the AIR volume for plotting
        pressures.push(P);
        times.push(time);
        thrusts.push(thrust);
        temperatures.push(temperature);

        // Update P_previous for next iteration
        P_previous = P;
        T_previous = temperature;
    }

    return {
        volumes,
        pressures,
        times,
        totalArea,
        thrusts,
        totalThrusts: totalThrust,
        temperatures
    };
}

function air(data) {
    // Unpack the arrays from the water phase
    const { temperatures, times, pressures, thrusts } = data;
    
    // Start from the final time, pressure, temperature from the water phase
    let timeCurrent     = times[times.length - 1];
    let pressureHigh = pressures[pressures.length - 1];
    let tempCurrent     = temperatures[temperatures.length - 1];
    
    // A small decrement in pressure each step
    // (smaller => more accurate but more iterations)
    const dP = 1000; 
    
    // We'll continue down to ambient (or until exactly atmospheric)
    while (pressureHigh > P_ATMOSPHERE) {
        // Next pressure step
        let pressureLow = pressureHigh - dP;
        if (pressureLow < P_ATMOSPHERE) {
            pressureLow = P_ATMOSPHERE;
        }
        
        // -- Compute next temperature (assuming adiabatic/isentropic) --
        //   For an isentropic process, T ~ P^((gamma-1)/gamma).
        //   We'll do a step-based approach:
        //        T_next = T_current * (P_next / P_current) ^ ((gamma-1)/gamma)
        let tempNext = tempCurrent * Math.pow(
            pressureLow / pressureHigh,
            (GAMMA - 1) / GAMMA
        );
        
        // --- Calculate densities and mass at current & next step ---
        const mAirCurrent = (pressureHigh / (GAS_CONSTANT * tempCurrent)) * V_BOTTLE;
        const mAirNext    = (pressureLow   / (GAS_CONSTANT * tempNext))    * V_BOTTLE;
        const deltaMass   = mAirCurrent - mAirNext; 
        // This is the small mass of air expelled going from P_current to P_next
        
        // --- Calculate an "average" temperature for the flow properties ---
        // (This is optional. We can also just use tempCurrent or tempNext.)
        let tempAvg = 0.5 * (tempCurrent + tempNext);
        
        // --- Calculate mass flow rate (approx for this step) ---
        const densityAvg   = pressureHigh / (GAS_CONSTANT * tempAvg);
        const veExp        = (GAMMA - 1) / GAMMA;
        const vePressure   = P_ATMOSPHERE / pressureHigh;
        const rhFormula    = 1 - Math.pow(vePressure, veExp);
        const lhFormula    = (2 * GAMMA * GAS_CONSTANT * tempAvg) / (GAMMA - 1);
        const velocityExhaust = Math.sqrt(lhFormula * rhFormula);
        
        const massFlowRate = densityAvg * NOZZLE_AREA * velocityExhaust;
        
        // --- Calculate thrust for this step ---
        // momentum thrust + pressure thrust
        const forceThrust =
            massFlowRate * velocityExhaust + 
            (pressureHigh - P_ATMOSPHERE) * NOZZLE_AREA;
        
        // --- Calculate the time for just this small mass to leave ---
        //   i.e., dt = deltaMass / massFlowRate
        //   (as opposed to massAir / massFlowRate, which would be "time to empty bottle")
        let dt = 0;
        if (massFlowRate > 1e-12 && deltaMass > 0) {
            dt = deltaMass / massFlowRate;
        } else {
            // If massFlowRate or deltaMass is extremely small, 
            // dt might blow up or be zero—handle gently:
            dt = 0; 
        }
        
        // Update the timeline
        timeCurrent += dt;
        console.log("pressureHigh:", pressureHigh, 
            "density", densityAvg,
            "exp", veExp,
            "ve pressure", vePressure,
            "pressureLow:", pressureLow,
            "temp:", tempAvg,
            "velocityExhaust:", velocityExhaust,
            "massFlowRate:", massFlowRate,
            "lh", lhFormula,
            "rh", rhFormula,
            "time", dt,
            "air mass", deltaMass
        );
        console.log("forceThrust:", forceThrust, "time", timeCurrent);

        // -- Push new data --
        pressures.push(pressureLow);
        temperatures.push(tempNext);
        thrusts.push(forceThrust);
        times.push(timeCurrent);
        
        // Advance to the next iteration
        pressureHigh = pressureLow;
        tempCurrent     = tempNext;
        
        // If we just reached atmospheric, break out
        if (pressureHigh <= P_ATMOSPHERE) {
            break;
        }
    }
    
    // Return updated arrays so plotAirThrustGraph() can access them
    return {
        thrusts,
        pressures,
        times,
        temperatures
    };
}




// Plot pressure vs volume graph
function plotPressureVolumeGraph(data) {
    const { volumes, pressures, totalArea } = data;

    const graphLayout = {
        autosize: true,
        margin: { l: 50, r: 20, t: 30, b: 40 },
        xaxis: {
            title: "Volume of Air (ml)",
            range: [0, 1500],
            tickformat: ",",
        },
        yaxis: { title: "Pressure (Pa)", range: [0, 500000] },
        annotations: [
            {
                xref: "paper",
                yref: "paper",
                x: 0.95,
                y: 0.05,
                xanchor: "right",
                yanchor: "bottom",
                text: `Total Area: ${totalArea.toFixed(2)} J`,
                showarrow: false,
                font: {
                    size: 12,
                    color: "black",
                },
            },
        ],
    };

    Plotly.newPlot("pressure-graph", [
        {
            x: volumes.map(v => v * 1000000), // Convert m³ to ml
            y: pressures,
            mode: "lines",
            name: "Pressure",
        },
        {
            x: [...volumes.map(v => v * 1000000), volumes[volumes.length - 1] * 1000000, volumes[0] * 1000000],
            y: [...pressures, 0, 0],
            fill: "tozeroy",
            type: "scatter",
            mode: "none",
            fillcolor: "rgba(0, 0, 255, 0.2)",
            name: "Area Under Curve",
        },
    ], { ...graphLayout, title: "Pressure vs Volume of Air" });
}

// Plot time vs volume of water graph
function plotTimeWaterVolumeGraph(data) {
    const { times } = data;
    const waterVolumes = data.volumes.map(v => V_BOTTLE - v); // Convert air volume to water volume

    const graphLayout = {
        autosize: true,
        margin: { l: 50, r: 20, t: 30, b: 40 },
        yaxis: {
            title: "Volume of Water (ml)",
            tickformat: ",",
            range: [0, 500],
        },
        xaxis: { title: "Time (s)", range: [0, 1] }, // Keep time on the y-axis
    };

    Plotly.newPlot("water-time-graph", [
        {
            y: waterVolumes.map(v => v * 1000000), // Convert m³ to ml
            x: times,
            mode: "lines",
            name: "Time",
        },
    ], { ...graphLayout, title: "Volume of Water Remaining vs Time" });
}

function plotTimeThrustGraph(data) {
    const { times, thrusts, totalThrusts } = data;

    const graphLayout = {
        autosize: true,
        margin: { l: 50, r: 20, t: 30, b: 40 },
        yaxis: {
            title: "Thrust (N)",
            tickformat: ",",
            range: [0, 60],
        },
        xaxis: { title: "Time (s)", range: [0, 1] }, // Keep time on the y-axis
        annotations: [
            {
                xref: "paper",
                yref: "paper",
                x: 0.95,
                y: 0.05,
                xanchor: "right",
                yanchor: "bottom",
                text: `Average Thrust: ${(totalThrusts / thrusts.length).toFixed(2)} N`,
                showarrow: false,
                font: {
                    size: 12,
                    color: "black",
                },
            },
        ],
    };

    Plotly.newPlot("thrust-time-graph", [
        {
            y: thrusts, // Convert m³ to ml
            x: times,
            mode: "lines",
            name: "Time",
        },
    ], { ...graphLayout, title: "Force Due to Thrust vs Time" });
}

function plotTempTimeGraph(data) {
    const { times, temperatures } = data;

    const graphLayout = {
        autosize: true,
        margin: { l: 50, r: 20, t: 30, b: 40 },
        yaxis: {
            title: "Temperature (K)",
            tickformat: ",",
            range: [0, 300],
        },
        xaxis: { title: "Time (s)", range: [0, 1] }, // Keep time on the y-axis
    };

    Plotly.newPlot("temperature-time-graph", [
        {
            y: temperatures, // Convert m³ to ml
            x: times,
            mode: "lines",
            name: "Time",
        },
    ], { ...graphLayout, title: "Temperature vs Time" });
}

function plotAirThrustGraph(data) {
    const { times, thrusts } = data;

    const graphLayout = {
        autosize: true,
        margin: { l: 50, r: 20, t: 30, b: 40 },
        yaxis: {
            title: "Thrust (N)",
            tickformat: ",",
            range: [0, 60],
        },
        xaxis: { title: "Time (s)", range: [0, 0.5] }, // Keep time on the y-axis
    };

    Plotly.newPlot("air-thrust-graph", [
        {
            y: thrusts, // Convert m³ to ml
            x: times,
            mode: "lines",
            name: "Time",
        },
    ], { ...graphLayout, title: "Thrust vs Time for Water and Air" });
}

// Initialize sliders and graphs
const psiSlider = document.getElementById("psi-slider");
const psiValue = document.getElementById("psi-value");
const waterVolumeSlider = document.getElementById("v-slider");
const waterVolumeValue = document.getElementById("v-value");

psiSlider.addEventListener("input", () => {
    psiValue.textContent = psiSlider.value;
    const data = simulate(parseFloat(psiSlider.value));
    plotPressureVolumeGraph(data); // Left graph remains unchanged
    plotTimeWaterVolumeGraph(data);
    plotTimeThrustGraph(data); // Updated right graph
    air(data);
    plotAirThrustGraph(data);
});

waterVolumeSlider.addEventListener("input", () => {
    waterVolumeValue.textContent = waterVolumeSlider.value;
    V_WATER = parseFloat(waterVolumeSlider.value) / 1000000; // Update water volume in m³
    const data = simulate(parseFloat(psiSlider.value));
    
    plotPressureVolumeGraph(data); // Left graph remains unchanged
    plotTimeWaterVolumeGraph(data); // Updated right graph
    plotTimeThrustGraph(data);
    plotTempTimeGraph(data);
    air(data);
    plotAirThrustGraph(data);
});

// Render initial graphs
const initialData = simulate(25);
plotPressureVolumeGraph(initialData); // Left graph remains unchanged
plotTimeWaterVolumeGraph(initialData); // Updated right graph
plotTimeThrustGraph(initialData);
plotTempTimeGraph(initialData);
air(initialData);
plotAirThrustGraph(initialData);

