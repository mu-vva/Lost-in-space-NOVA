# Lost-in-space-NOVA

Case adaptive satellite imaging scheduler developed for the 418 AEON Hackathon

# Satellite Imaging Planner

### Team Nova — 418 AEON TakeMe2Space Hackathon

**Team Members**
- Muvva Sirimalli
- Satya Pranavi Vemula
- Kundana Priya

---

## Overview

This project implements a **case-adaptive imaging scheduler** for an Earth observation satellite.
The objective is to maximize **Area of Interest (AOI) coverage** while respecting spacecraft constraints and improving control-effort efficiency.
The solution uses different scheduling strategies depending on the visibility and geometric constraints of each scenario.

---

## Key Idea

The scoring function prioritizes:

1. **Coverage (C)** as the dominant objective
2. **Control-effort efficiency (η_E)** as a secondary objective

A key observation was that high peak angular velocity during fast slews can result in a large momentum change (ΔH_used), reducing control-effort efficiency.

Our strategy therefore focuses on:

> **Spreading imaging operations over time to reduce peak angular velocity while maintaining high coverage.**

---

## Approach

### Case-Adaptive Planning

Different scheduling strategies are used depending on the available imaging window.

#### Case 1 & Case 2 — Broad Visibility

For scenarios with broader visibility windows, the planner uses **time-slot scheduling**.

- Frames are distributed across the available pass.
- A larger `min_gap` (approximately 2.5 seconds) is used.
- This reduces the need for rapid consecutive slews.
- Lower peak slew rates improve control-effort efficiency.

#### Case 3 — Narrow Visibility

For the more constrained visibility window, the planner uses a **greedy best-time scheduling strategy**.

- A smaller `min_gap` (approximately 0.9 seconds) is used.
- The scheduler prioritizes capturing as many frames as possible within the limited window.
- Coverage is prioritized over control-effort efficiency when the geometry makes spreading the observations impractical.

---

## Path Optimization

For Cases 1 and 2, path optimization is performed before scheduling.

The pipeline consists of:

1. **Strip-snake ordering**
   - Provides an initial spatially coherent ordering of imaging targets.

2. **Nearest-neighbour reordering**
   - Reduces angular jumps between consecutive targets.

3. **Local path improvement**
   - Further reduces unnecessary movement along the imaging sequence.

This reduces total rotation and helps lower the control effort required by the spacecraft.

---

## Attitude Strategy

Each imaging frame follows a fixed attitude sequence:

| Stage | Duration |
|---|---:|
| Pre-hold | 50 ms |
| Shutter | 120 ms |
| Post-hold | 50 ms |

During the shutter interval, the spacecraft maintains a constant attitude to satisfy the smear constraint.

Between imaging frames, the attitude transitions using smooth interpolation.

---

## Constraint Handling

The planner enforces the key spacecraft constraints throughout the schedule.

### Smear Constraint

The angular velocity constraint is:

`|ω| ≤ 0.05°/s`

The spacecraft maintains a constant attitude during image capture to avoid image smear.

### Off-Nadir Constraint

The planner maintains:

`Off-nadir ≤ 60°`

with a safety margin.

### Shutter Overlap

Imaging windows are prevented from overlapping through the use of the `min_gap` parameter.

### Monotonic Timeline

The generated attitude timeline remains chronologically ordered throughout the schedule.

### Wheel Limits

Reaction-wheel constraints are addressed indirectly by:

- limiting slew rates
- spreading maneuvers over time
- avoiding unnecessarily rapid attitude changes

---

## Key Insight

The main insight from the scheduling experiments was:

> **η_E is improved not simply by minimizing motion, but by reducing peak angular velocity through appropriate time spacing.**

A schedule with slightly more total movement can therefore achieve better control-effort efficiency if the movement is distributed over a longer period rather than concentrated into rapid slews.

---

## Expected Performance

The solution produced the following expected behavior across the test cases:

| Scenario | Approx. Frames | Expected η_E | Strategy |
|---|---:|---:|---|
| Case 1 | ~49 | ~0.3–0.4 | Time-slot scheduling |
| Case 2 | ~49 | ~0.3–0.4 | Time-slot scheduling |
| Case 3 | ~45–47 | ~0 | Greedy best-time scheduling |

Case 3 prioritizes coverage because the narrow visibility window limits the ability to distribute maneuvers over time.

---

## Overall Strategy

The solution combines **greedy scheduling, spatial path optimization, and time-distributed planning**.

The overall objective is to balance:

- **High AOI coverage** — primary objective
- **Controlled spacecraft motion** — secondary objective
- **Constraint compliance** — required throughout the schedule

The planner therefore adapts its scheduling strategy to the geometry and available visibility of each scenario rather than using one fixed strategy for all cases.

---

## Repository Contents

```text
satellite-imaging-planner/
│
├── README.md
└── solution.[file-extension]
