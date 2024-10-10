/* Copyright 2023 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */

#include "HybridPICModel.H"

#include "EmbeddedBoundary/Enabled.H"
#include "Fields.H"
#include "WarpX.H"
#include "Utils/Algorithms/LinearInterpolation.H"

using namespace amrex;
using warpx::fields::FieldType;


HybridPICModel::HybridPICModel ()
{
    ReadParameters();
}

void HybridPICModel::ReadParameters ()
{
    const ParmParse pp_hybrid("hybrid_pic_model");

    // The B-field update is subcycled to improve stability - the number
    // of sub steps can be specified by the user (defaults to 50).
    utils::parser::queryWithParser(pp_hybrid, "substeps", m_substeps);

    // The hybrid model requires an electron temperature, reference density
    // and exponent to be given. These values will be used to calculate the
    // electron pressure according to p = n0 * Te * (n/n0)^gamma
    utils::parser::queryWithParser(pp_hybrid, "gamma", m_gamma);
    if (!pp_hybrid.query("elec_temp", m_elec_temp_expression)) {
        Abort("hybrid_pic_model.elec_temp must be specified when using the hybrid solver");
    }
    if (!pp_hybrid.query("electron_temperature_init_style", m_elec_temp_style)) {
        Abort("hybrid_pic_model.electron_temperature_init_style must be specified "
              "when using the hybrid solver");
    }

    if (m_elec_temp_style == "read_from_file") {
        pp_hybrid.get("read_Te_field_from_path", m_elec_temp_field_path);
    }

    const bool n0_ref_given = utils::parser::queryWithParser(pp_hybrid, "n0_ref", m_n0_ref);
    if (m_gamma != 1.0 && !n0_ref_given) {
        Abort("hybrid_pic_model.n0_ref should be specified if hybrid_pic_model.gamma != 1");
    }

    pp_hybrid.query("plasma_resistivity(rho,J)", m_eta_expression);
    utils::parser::queryWithParser(pp_hybrid, "n_floor", m_n_floor);

    utils::parser::queryWithParser(pp_hybrid, "plasma_hyper_resistivity", m_eta_h);

    // external currents
    pp_hybrid.query("Jx_external_grid_function(x,y,z,t)", m_Jx_ext_grid_function);
    pp_hybrid.query("Jy_external_grid_function(x,y,z,t)", m_Jy_ext_grid_function);
    pp_hybrid.query("Jz_external_grid_function(x,y,z,t)", m_Jz_ext_grid_function);

    pp_hybrid.query("J_external_init_style", m_J_ext_grid_style);

    if (m_J_ext_grid_style == "read_from_file") {
        pp_hybrid.get("read_j_fields_from_path", m_external_j_fields_path);
    }

    // external magnetic field
    pp_hybrid.query("Bx_external_grid_function(x,y,z,t)", m_Bx_ext_grid_function);
    pp_hybrid.query("By_external_grid_function(x,y,z,t)", m_By_ext_grid_function);
    pp_hybrid.query("Bz_external_grid_function(x,y,z,t)", m_Bz_ext_grid_function);

    pp_hybrid.query("B_external_init_style", m_B_ext_grid_style);

    if (m_B_ext_grid_style == "read_from_file"){
        pp_hybrid.get("read_b_fields_from_path", m_external_b_fields_path);
    }
}

void HybridPICModel::AllocateLevelMFs (ablastr::fields::MultiFabRegister & fields,
                                       int lev, const BoxArray& ba, const DistributionMapping& dm,
                                       const int ncomps, const IntVect& ngEB, const IntVect& ngJ, const IntVect& ngRho,
                                       const IntVect& Bx_nodal_flag,
                                       const IntVect& By_nodal_flag,
                                       const IntVect& Bz_nodal_flag,
                                       const IntVect& jx_nodal_flag,
                                       const IntVect& jy_nodal_flag,
                                       const IntVect& jz_nodal_flag,
                                       const IntVect& rho_nodal_flag)
{
    using ablastr::fields::Direction;

    // The "hybrid_electron_pressure_fp" multifab stores the electron pressure calculated
    // from the specified equation of state.
    fields.alloc_init(FieldType::hybrid_electron_pressure_fp,
        lev, amrex::convert(ba, rho_nodal_flag),
        dm, ncomps, ngRho, 0.0_rt);

    // The "hybrid_electron_temperature_fp" multifab stores the electron termperature
    fields.alloc_init(FieldType::hybrid_electron_temperature_fp,
        lev, amrex::convert(ba, rho_nodal_flag),
        dm, ncomps, ngRho, 0.0_rt);

    // The "hybrid_rho_fp_temp" multifab is used to store the ion charge density
    // interpolated or extrapolated to appropriate timesteps.
    fields.alloc_init(FieldType::hybrid_rho_fp_temp,
        lev, amrex::convert(ba, rho_nodal_flag),
        dm, ncomps, ngRho, 0.0_rt);

    // The "hybrid_current_fp_temp" multifab is used to store the ion current density
    // interpolated or extrapolated to appropriate timesteps.
    fields.alloc_init(FieldType::hybrid_current_fp_temp, Direction{0},
        lev, amrex::convert(ba, jx_nodal_flag),
        dm, ncomps, ngJ, 0.0_rt);
    fields.alloc_init(FieldType::hybrid_current_fp_temp, Direction{1},
        lev, amrex::convert(ba, jy_nodal_flag),
        dm, ncomps, ngJ, 0.0_rt);
    fields.alloc_init(FieldType::hybrid_current_fp_temp, Direction{2},
        lev, amrex::convert(ba, jz_nodal_flag),
        dm, ncomps, ngJ, 0.0_rt);

    // The "hybrid_current_fp_plasma" multifab stores the total plasma current calculated
    // as the curl of B minus any external current.
    fields.alloc_init(FieldType::hybrid_current_fp_plasma, Direction{0},
        lev, amrex::convert(ba, jx_nodal_flag),
        dm, ncomps, ngJ, 0.0_rt);
    fields.alloc_init(FieldType::hybrid_current_fp_plasma, Direction{1},
        lev, amrex::convert(ba, jy_nodal_flag),
        dm, ncomps, ngJ, 0.0_rt);
    fields.alloc_init(FieldType::hybrid_current_fp_plasma, Direction{2},
        lev, amrex::convert(ba, jz_nodal_flag),
        dm, ncomps, ngJ, 0.0_rt);

    // the external current density multifab matches the current staggering and
    // one ghost cell is used since we interpolate the current to a nodal grid
    fields.alloc_init(FieldType::hybrid_current_fp_external, Direction{0},
        lev, amrex::convert(ba, jx_nodal_flag),
        dm, ncomps, IntVect(1), 0.0_rt);
    fields.alloc_init(FieldType::hybrid_current_fp_external, Direction{1},
        lev, amrex::convert(ba, jy_nodal_flag),
        dm, ncomps, IntVect(1), 0.0_rt);
    fields.alloc_init(FieldType::hybrid_current_fp_external, Direction{2},
        lev, amrex::convert(ba, jz_nodal_flag),
        dm, ncomps, IntVect(1), 0.0_rt);

    fields.alloc_init(FieldType::hybrid_bfield_fp_external, Direction{0},
        lev, amrex::convert(ba, Bx_nodal_flag), dm, ncomps, ngEB, 0.0_rt);
    fields.alloc_init(FieldType::hybrid_bfield_fp_external, Direction{1},
        lev, amrex::convert(ba, By_nodal_flag), dm, ncomps, ngEB, 0.0_rt);
    fields.alloc_init(FieldType::hybrid_bfield_fp_external, Direction{2},
        lev, amrex::convert(ba, Bz_nodal_flag), dm, ncomps, ngEB, 0.0_rt);

#ifdef WARPX_DIM_RZ
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        (ncomps == 1),
        "Ohm's law solver only support m = 0 azimuthal mode at present.");
#endif
}

void HybridPICModel::InitData ()
{
    m_elec_temp_parser = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(m_elec_temp_expression, {"x","y","z"}));
    m_elec_temp = m_elec_temp_parser->compile<3>();

    m_resistivity_parser = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(m_eta_expression, {"rho","J"}));
    m_eta = m_resistivity_parser->compile<2>();
    const std::set<std::string> resistivity_symbols = m_resistivity_parser->symbols();
    m_resistivity_has_J_dependence += resistivity_symbols.count("J");

    m_J_external_parser[0] = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(m_Jx_ext_grid_function,{"x","y","z","t"}));
    m_J_external_parser[1] = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(m_Jy_ext_grid_function,{"x","y","z","t"}));
    m_J_external_parser[2] = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(m_Jz_ext_grid_function,{"x","y","z","t"}));
    m_J_external[0] = m_J_external_parser[0]->compile<4>();
    m_J_external[1] = m_J_external_parser[1]->compile<4>();
    m_J_external[2] = m_J_external_parser[2]->compile<4>();

    // check if the external current parsers depend on time
    for (int i=0; i<3; i++) {
        const std::set<std::string> J_ext_symbols = m_J_external_parser[i]->symbols();
        m_external_jfield_has_time_dependence += J_ext_symbols.count("t");
    }

    m_B_external_parser[0] = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(m_Bx_ext_grid_function,{"x","y","z","t"}));
    m_B_external_parser[1] = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(m_By_ext_grid_function,{"x","y","z","t"}));
    m_B_external_parser[2] = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(m_Bz_ext_grid_function,{"x","y","z","t"}));
    m_B_external[0] = m_B_external_parser[0]->compile<4>();
    m_B_external[1] = m_B_external_parser[1]->compile<4>();
    m_B_external[2] = m_B_external_parser[2]->compile<4>();

    // check if the external b field parsers depend on time
    for (int i=0; i<3; i++) {
        const std::set<std::string> B_ext_symbols = m_B_external_parser[i]->symbols();
        m_external_bfield_has_time_dependence += B_ext_symbols.count("t");
    }

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !m_external_bfield_has_time_dependence,
        "The hybrid-PIC algorithm does not work with time dependent external magnetic fields."
    );

    auto & warpx = WarpX::GetInstance();
    using ablastr::fields::Direction;

    // Get the grid staggering of the fields involved in calculating E
    amrex::IntVect Jx_stag = warpx.m_fields.get(FieldType::current_fp, Direction{0}, 0)->ixType().toIntVect();
    amrex::IntVect Jy_stag = warpx.m_fields.get(FieldType::current_fp, Direction{1}, 0)->ixType().toIntVect();
    amrex::IntVect Jz_stag = warpx.m_fields.get(FieldType::current_fp, Direction{2}, 0)->ixType().toIntVect();
    amrex::IntVect Bx_stag = warpx.m_fields.get(FieldType::Bfield_fp, Direction{0}, 0)->ixType().toIntVect();
    amrex::IntVect By_stag = warpx.m_fields.get(FieldType::Bfield_fp, Direction{1}, 0)->ixType().toIntVect();
    amrex::IntVect Bz_stag = warpx.m_fields.get(FieldType::Bfield_fp, Direction{2}, 0)->ixType().toIntVect();
    amrex::IntVect Ex_stag = warpx.m_fields.get(FieldType::Efield_fp, Direction{0}, 0)->ixType().toIntVect();
    amrex::IntVect Ey_stag = warpx.m_fields.get(FieldType::Efield_fp, Direction{1}, 0)->ixType().toIntVect();
    amrex::IntVect Ez_stag = warpx.m_fields.get(FieldType::Efield_fp, Direction{2}, 0)->ixType().toIntVect();

    // Check that the grid types are appropriate
    const bool appropriate_grids = (
#if   defined(WARPX_DIM_1D_Z)
        // AMReX convention: x = missing dimension, y = missing dimension, z = only dimension
        Ex_stag == IntVect(1) && Ey_stag == IntVect(1) && Ez_stag == IntVect(0) &&
        Bx_stag == IntVect(0) && By_stag == IntVect(0) && Bz_stag == IntVect(1) &&
#elif   defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        // AMReX convention: x = first dimension, y = missing dimension, z = second dimension
        Ex_stag == IntVect(0,1) && Ey_stag == IntVect(1,1) && Ez_stag == IntVect(1,0) &&
        Bx_stag == IntVect(1,0) && By_stag == IntVect(0,0) && Bz_stag == IntVect(0,1) &&
#elif defined(WARPX_DIM_3D)
        Ex_stag == IntVect(0,1,1) && Ey_stag == IntVect(1,0,1) && Ez_stag == IntVect(1,1,0) &&
        Bx_stag == IntVect(1,0,0) && By_stag == IntVect(0,1,0) && Bz_stag == IntVect(0,0,1) &&
#endif
        Jx_stag == Ex_stag && Jy_stag == Ey_stag && Jz_stag == Ez_stag
    );
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        appropriate_grids,
        "Ohm's law E-solve only works with staggered (Yee) grids.");

    // copy data to device
    for ( int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Jx_IndexType[idim]    = Jx_stag[idim];
        Jy_IndexType[idim]    = Jy_stag[idim];
        Jz_IndexType[idim]    = Jz_stag[idim];
        Bx_IndexType[idim]    = Bx_stag[idim];
        By_IndexType[idim]    = By_stag[idim];
        Bz_IndexType[idim]    = Bz_stag[idim];
        Ex_IndexType[idim]    = Ex_stag[idim];
        Ey_IndexType[idim]    = Ey_stag[idim];
        Ez_IndexType[idim]    = Ez_stag[idim];
    }

    // Below we set all the unused dimensions to have nodal values for J, B & E
    // since these values will be interpolated onto a nodal grid - if this is
    // not done the Interp function returns nonsense values.
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ) || defined(WARPX_DIM_1D_Z)
    Jx_IndexType[2]    = 1;
    Jy_IndexType[2]    = 1;
    Jz_IndexType[2]    = 1;
    Bx_IndexType[2]    = 1;
    By_IndexType[2]    = 1;
    Bz_IndexType[2]    = 1;
    Ex_IndexType[2]    = 1;
    Ey_IndexType[2]    = 1;
    Ez_IndexType[2]    = 1;
#endif
#if defined(WARPX_DIM_1D_Z)
    Jx_IndexType[1]    = 1;
    Jy_IndexType[1]    = 1;
    Jz_IndexType[1]    = 1;
    Bx_IndexType[1]    = 1;
    By_IndexType[1]    = 1;
    Bz_IndexType[1]    = 1;
    Ex_IndexType[1]    = 1;
    Ey_IndexType[1]    = 1;
    Ez_IndexType[1]    = 1;
#endif

    // Initialize external current - note that this approach skips the check
    // if the current is time dependent which is what needs to be done to
    // write time independent fields on the first step.
    for (int lev = 0; lev <= warpx.finestLevel(); ++lev) {
        ablastr::fields::ScalarField electron_temperature_fp = warpx.m_fields.get(FieldType::hybrid_electron_temperature_fp, lev);
        if (m_elec_temp_style == "parse_Te_ext_grid_function") {
            ComputeExternalScalarFieldOnGridUsingParser(
                FieldType::hybrid_electron_temperature_fp,
                m_elec_temp,
                lev, PatchType::fine
            );
        } else if (m_elec_temp_style == "read_from_file") {
            ReadScalarFieldFromFile(
                m_elec_temp_field_path,
                *electron_temperature_fp);
        }

        if (m_J_ext_grid_style == "parse_j_ext_grid_function") {
            warpx.ComputeExternalFieldOnGridUsingParser(
                FieldType::hybrid_current_fp_external,
                m_J_external[0],
                m_J_external[1],
                m_J_external[2],
                lev, PatchType::fine, 'e',
                warpx.m_fields.get_alldirs(FieldType::edge_lengths, lev),
                warpx.m_fields.get_alldirs(FieldType::face_areas, lev));
        } else if (m_J_ext_grid_style == "read_from_file") {
            warpx.ReadExternalFieldFromFile(
                m_external_j_fields_path,
                warpx.m_fields.get(FieldType::hybrid_current_fp_external, Direction{0}, lev),
                "J", "x");
            warpx.ReadExternalFieldFromFile(
                m_external_j_fields_path,
                warpx.m_fields.get(FieldType::hybrid_current_fp_external, Direction{1}, lev),
                "J", "y");
            warpx.ReadExternalFieldFromFile(
                m_external_j_fields_path,
                warpx.m_fields.get(FieldType::hybrid_current_fp_external, Direction{2}, lev),
                "J", "z");
        }

        if (m_B_ext_grid_style == "parse_b_ext_grid_function") {
            warpx.ComputeExternalFieldOnGridUsingParser(
                FieldType::hybrid_bfield_fp_external,
                m_B_external[0],
                m_B_external[1],
                m_B_external[2],
                lev, PatchType::fine, 'f',
                warpx.m_fields.get_alldirs(FieldType::edge_lengths, lev),
                warpx.m_fields.get_alldirs(FieldType::face_areas, lev));
        } else if (m_B_ext_grid_style == "read_from_file") {
            warpx.ReadExternalFieldFromFile(
                m_external_b_fields_path,
                warpx.m_fields.get(FieldType::hybrid_bfield_fp_external, Direction{0}, lev),
                "B", "x");
            warpx.ReadExternalFieldFromFile(
                m_external_b_fields_path,
                warpx.m_fields.get(FieldType::hybrid_bfield_fp_external, Direction{1}, lev),
                "B", "y");
            warpx.ReadExternalFieldFromFile(
                m_external_b_fields_path,
                warpx.m_fields.get(FieldType::hybrid_bfield_fp_external, Direction{2}, lev),
                "B", "z");
        }
    }
}

void HybridPICModel::GetCurrentExternal ()
{
    if (!m_external_jfield_has_time_dependence) { return; }

    auto& warpx = WarpX::GetInstance();
    for (int lev = 0; lev <= warpx.finestLevel(); ++lev)
    {
        warpx.ComputeExternalFieldOnGridUsingParser(
            FieldType::hybrid_current_fp_external,
            m_J_external[0],
            m_J_external[1],
            m_J_external[2],
            lev, PatchType::fine, 'e',
            warpx.m_fields.get_alldirs(FieldType::edge_lengths, lev),
            warpx.m_fields.get_alldirs(FieldType::face_areas, lev));
    }
}

void HybridPICModel::CalculatePlasmaCurrent (
    ablastr::fields::MultiLevelVectorField const& Bfield,
    ablastr::fields::MultiLevelVectorField const& edge_lengths)
{
    auto& warpx = WarpX::GetInstance();
    for (int lev = 0; lev <= warpx.finestLevel(); ++lev)
    {
        CalculatePlasmaCurrent(Bfield[lev], edge_lengths[lev], lev);
    }
}

void HybridPICModel::CalculatePlasmaCurrent (
    ablastr::fields::VectorField const& Bfield,
    ablastr::fields::VectorField const& edge_lengths,
    const int lev)
{
    WARPX_PROFILE("HybridPICModel::CalculatePlasmaCurrent()");

    auto& warpx = WarpX::GetInstance();
    ablastr::fields::VectorField current_fp_plasma = warpx.m_fields.get_alldirs(FieldType::hybrid_current_fp_plasma, lev);
    warpx.get_pointer_fdtd_solver_fp(lev)->CalculateCurrentAmpere(
        current_fp_plasma, Bfield, edge_lengths, lev
    );

    // we shouldn't apply the boundary condition to J since J = J_i - J_e but
    // the boundary correction was already applied to J_i and the B-field
    // boundary ensures that J itself complies with the boundary conditions, right?
    // ApplyJfieldBoundary(lev, Jfield[0].get(), Jfield[1].get(), Jfield[2].get());
    for (int i=0; i<3; i++) { current_fp_plasma[i]->FillBoundary(warpx.Geom(lev).periodicity()); }

    // Subtract external current from "Ampere" current calculated above. Note
    // we need to include 1 ghost cell since later we will interpolate the
    // plasma current to a nodal grid.
    ablastr::fields::VectorField current_fp_external = warpx.m_fields.get_alldirs(FieldType::hybrid_current_fp_external, lev);
    for (int i=0; i<3; i++) {
        current_fp_plasma[i]->minus(*current_fp_external[i], 0, 1, 1);
    }

}
void HybridPICModel::ComputeExternalScalarFieldOnGridUsingParser (
    warpx::fields::FieldType field,
    amrex::ParserExecutor<3> const& scalar_parser,
    int lev, PatchType patch_type)
{

    auto& warpx = WarpX::GetInstance();
    const amrex::Geometry& geom = warpx.Geom(lev);
    auto dx_lev = geom.CellSizeArray();
    const RealBox& real_box = geom.ProbDomain();

    amrex::IntVect refratio = (lev > 0 ) ? warpx.RefRatio(lev-1) : amrex::IntVect(1);
    if (patch_type == PatchType::coarse) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            dx_lev[idim] = dx_lev[idim] * refratio[idim];
        }
    }

    ablastr::fields::ScalarField mf = warpx.m_fields.get(field, lev);

    const amrex::IntVect nodal_flag = mf->ixType().toIntVect();

    // Loop over boxes
    for (amrex::MFIter mfi(*mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const& scalar_arr = mf->array(mfi);

        // Start ParallelFor
        amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                // Shift required in the x-, y-, or z- position
                // depending on the index type of the multifab
#if defined(WARPX_DIM_1D_Z)
                const amrex::Real x = 0._rt;
                const amrex::Real y = 0._rt;
                const amrex::Real fac_z = (1._rt - nodal_flag[0]) * dx_lev[0] * 0.5_rt;
                const amrex::Real z = i*dx_lev[0] + real_box.lo(0) + fac_z;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                const amrex::Real fac_x = (1._rt - nodal_flag[0]) * dx_lev[0] * 0.5_rt;
                const amrex::Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
                const amrex::Real y = 0._rt;
                const amrex::Real fac_z = (1._rt - nodal_flag[1]) * dx_lev[1] * 0.5_rt;
                const amrex::Real z = j*dx_lev[1] + real_box.lo(1) + fac_z;
#else
                const amrex::Real fac_x = (1._rt - nodal_flag[0]) * dx_lev[0] * 0.5_rt;
                const amrex::Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
                const amrex::Real fac_y = (1._rt - nodal_flag[1]) * dx_lev[1] * 0.5_rt;
                const amrex::Real y = j*dx_lev[1] + real_box.lo(1) + fac_y;
                const amrex::Real fac_z = (1._rt - nodal_flag[2]) * dx_lev[2] * 0.5_rt;
                const amrex::Real z = k*dx_lev[2] + real_box.lo(2) + fac_z;
#endif
                // Assign the value to the scalar field array
                scalar_arr(i, j, k) = scalar_parser(x, y, z) * PhysConst::q_e;
            });
    }
}

void HybridPICModel::ReadScalarFieldFromFile (
    const std::string& scalar_file,
    amrex::MultiFab& scalar_field)
{
    // Open the file
    std::ifstream infile(scalar_file);
    if (!infile.is_open()) {
        amrex::Abort("Failed to open file: " + scalar_file);
    }

    auto& warpx = WarpX::GetInstance();
    const amrex::Geometry& geom = warpx.Geom(0);

    // Read the data into a host vector
    const auto& domain = geom.Domain();
    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();
    const int nx = domain.length(0);
    const int ny = domain.length(1);
    const int nz = domain.length(2);
    const size_t total_cells = nx * ny * nz;

    std::vector<double> scalar_data_host(total_cells);
    for (size_t i = 0; i < total_cells; ++i) {
        infile >> scalar_data_host[i];
    }
    infile.close();

    // Copy data to GPU
    amrex::Gpu::DeviceVector<double> scalar_data_gpu(total_cells);
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, scalar_data_host.begin(), scalar_data_host.end(), scalar_data_gpu.begin());

    // Get a raw pointer to the GPU data
    double* scalar_data_ptr = scalar_data_gpu.data();

    // Loop over boxes
    for (amrex::MFIter mfi(scalar_field, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const& temp_arr = scalar_field.array(mfi);

        // Start ParallelFor
        amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                // Compute the physical coordinates of the grid point
                amrex::Real x0, x1, x2;
                if (bx.type(0) == amrex::IndexType::CellIndex::NODE) {
                    x0 = prob_lo[0] + i * dx[0];
                } else {
                    x0 = prob_lo[0] + i * dx[0] + 0.5_rt * dx[0];
                }
                if (bx.type(1) == amrex::IndexType::CellIndex::NODE) {
                    x1 = prob_lo[1] + j * dx[1];
                } else {
                    x1 = prob_lo[1] + j * dx[1] + 0.5_rt * dx[1];
                }
#if defined(WARPX_DIM_3D)
                if (bx.type(2) == amrex::IndexType::CellIndex::NODE) {
                    x2 = prob_lo[2] + k * dx[2];
                } else {
                    x2 = prob_lo[2] + k * dx[2] + 0.5_rt * dx[2];
                }
#endif

                // Map the physical coordinates to the OpenPMD grid
                int ix = static_cast<int>((x0 - prob_lo[0]) / dx[0]);
                int iy = static_cast<int>((x1 - prob_lo[1]) / dx[1]);
#if defined(WARPX_DIM_3D)
                int iz = static_cast<int>((x2 - prob_lo[2]) / dx[2]);
#else
                int iz = k; // In 2D, k is just 0
#endif

                // Get coordinates of external grid point
                amrex::Real xx0 = prob_lo[0] + ix * dx[0];
                amrex::Real xx1 = prob_lo[1] + iy * dx[1];
#if defined(WARPX_DIM_3D)
                amrex::Real xx2 = prob_lo[2] + iz * dx[2];
#endif

                // Interpolation
#if defined(WARPX_DIM_RZ)
                const amrex::Array4<double> fc_array(scalar_data_ptr, {0, 0, 0}, {nx, ny, nz}, 1);
                const double
                    f00 = fc_array(0, iz, ix),
                    f01 = fc_array(0, iz, ix + 1),
                    f10 = fc_array(0, iz + 1, ix),
                    f11 = fc_array(0, iz + 1, ix + 1);
                temp_arr(i, j, k) = static_cast<amrex::Real>(utils::algorithms::bilinear_interp<double>
                    (xx0, xx0 + dx[0], xx1, xx1 + dx[1],
                     f00, f01, f10, f11,
                     x0, x1));
#elif defined(WARPX_DIM_3D)
                const amrex::Array4<double> fc_array(scalar_data_ptr, {0, 0, 0}, {nx, ny, nz}, 1);
                const double
                    f000 = fc_array(iz, iy, ix),
                    f001 = fc_array(iz, iy, ix + 1),
                    f010 = fc_array(iz, iy + 1, ix),
                    f011 = fc_array(iz, iy + 1, ix + 1),
                    f100 = fc_array(iz + 1, iy, ix),
                    f101 = fc_array(iz + 1, iy, ix + 1),
                    f110 = fc_array(iz + 1, iy + 1, ix),
                    f111 = fc_array(iz + 1, iy + 1, ix + 1);
                temp_arr(i, j, k) = static_cast<amrex::Real>(utils::algorithms::trilinear_interp<double>
                    (xx0, xx0 + dx[0], xx1, xx1 + dx[1], xx2, xx2 + dx[2],
                     f000, f001, f010, f011, f100, f101, f110, f111,
                     x0, x1, x2));
#endif
                temp_arr(i, j, k) = temp_arr(i, j, k) * PhysConst::q_e;
            });
    }
}

void HybridPICModel::HybridPICSolveE (
    ablastr::fields::MultiLevelVectorField const& Efield,
    ablastr::fields::MultiLevelVectorField const& Jfield,
    ablastr::fields::MultiLevelVectorField const& Bfield,
    ablastr::fields::MultiLevelScalarField const& rhofield,
    ablastr::fields::MultiLevelVectorField const& edge_lengths,
    const bool solve_for_Faraday) const
{
    auto& warpx = WarpX::GetInstance();
    for (int lev = 0; lev <= warpx.finestLevel(); ++lev)
    {
        HybridPICSolveE(
            Efield[lev], Jfield[lev], Bfield[lev], *rhofield[lev],
            edge_lengths[lev], lev, solve_for_Faraday
        );
    }
}

void HybridPICModel::HybridPICSolveE (
    ablastr::fields::VectorField const& Efield,
    ablastr::fields::VectorField const& Jfield,
    ablastr::fields::VectorField const& Bfield,
    amrex::MultiFab const& rhofield,
    ablastr::fields::VectorField const& edge_lengths,
    const int lev, const bool solve_for_Faraday) const
{
    WARPX_PROFILE("WarpX::HybridPICSolveE()");

    HybridPICSolveE(
        Efield, Jfield, Bfield, rhofield, edge_lengths, lev,
        PatchType::fine, solve_for_Faraday
    );
    if (lev > 0)
    {
        amrex::Abort(Utils::TextMsg::Err(
        "HybridPICSolveE: Only one level implemented for hybrid-PIC solver."));
    }
}

void HybridPICModel::HybridPICSolveE (
    ablastr::fields::VectorField const& Efield,
    ablastr::fields::VectorField const& Jfield,
    ablastr::fields::VectorField const& Bfield,
    amrex::MultiFab const& rhofield,
    ablastr::fields::VectorField const& edge_lengths,
    const int lev, PatchType patch_type,
    const bool solve_for_Faraday) const
{
    auto& warpx = WarpX::GetInstance();

    ablastr::fields::VectorField current_fp_plasma = warpx.m_fields.get_alldirs(FieldType::hybrid_current_fp_plasma, lev);
    const ablastr::fields::ScalarField electron_pressure_fp = warpx.m_fields.get(FieldType::hybrid_electron_pressure_fp, lev);

    // Solve E field in regular cells
    warpx.get_pointer_fdtd_solver_fp(lev)->HybridPICSolveE(
        Efield, current_fp_plasma, Jfield, Bfield, rhofield,
        *electron_pressure_fp, edge_lengths, lev, this, solve_for_Faraday
    );
    warpx.ApplyEfieldBoundary(lev, patch_type);
}

void HybridPICModel::CalculateElectronPressure() const
{
    auto& warpx = WarpX::GetInstance();
    for (int lev = 0; lev <= warpx.finestLevel(); ++lev)
    {
        CalculateElectronPressure(lev);
    }
}

void HybridPICModel::CalculateElectronPressure(const int lev) const
{
    WARPX_PROFILE("WarpX::CalculateElectronPressure()");

    auto& warpx = WarpX::GetInstance();
    ablastr::fields::ScalarField electron_pressure_fp = warpx.m_fields.get(FieldType::hybrid_electron_pressure_fp, lev);
    ablastr::fields::ScalarField electron_temperature_fp = warpx.m_fields.get(FieldType::hybrid_electron_temperature_fp, lev);
    ablastr::fields::ScalarField rho_fp = warpx.m_fields.get(FieldType::rho_fp, lev);

    // Calculate the electron pressure using rho^{n+1}.

    FillElectronPressureMF(
        *electron_pressure_fp,
        *electron_temperature_fp,
        *rho_fp
    );


    warpx.ApplyElectronPressureBoundary(lev, PatchType::fine);
    electron_pressure_fp->FillBoundary(warpx.Geom(lev).periodicity());
}

void HybridPICModel::FillElectronPressureMF (
    amrex::MultiFab& Pe_field,
    amrex::MultiFab const& Te_field,
    amrex::MultiFab const& rho_field
) const
{
    const auto n0_ref = m_n0_ref;
    const auto gamma = m_gamma;

    // Loop through the grids, and over the tiles within each grid
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(Pe_field, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Extract field data for this grid/tile
        Array4<Real const> const& rho = rho_field.const_array(mfi);
        Array4<Real const> const& Te = Te_field.const_array(mfi);
        Array4<Real> const& Pe = Pe_field.array(mfi);

        // Extract tileboxes for which to loop
        const Box& tilebox  = mfi.tilebox();

        ParallelFor(tilebox, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Pe(i, j, k) = ElectronPressure::get_pressure(
                n0_ref, Te(i, j, k), gamma, rho(i, j, k)
            );
        });
    }
}

void HybridPICModel::BfieldEvolveRK (
    ablastr::fields::MultiLevelVectorField const& Bfield,
    ablastr::fields::MultiLevelVectorField const& Efield,
    ablastr::fields::MultiLevelVectorField const& Jfield,
    ablastr::fields::MultiLevelScalarField const& rhofield,
    ablastr::fields::MultiLevelVectorField  const& edge_lengths,
    amrex::Real dt, DtType dt_type,
    IntVect ng, std::optional<bool> nodal_sync )
{
    auto& warpx = WarpX::GetInstance();
    for (int lev = 0; lev <= warpx.finestLevel(); ++lev)
    {
        BfieldEvolveRK(
            Bfield, Efield, Jfield, rhofield, edge_lengths, dt, lev, dt_type,
            ng, nodal_sync
        );
    }
}

void HybridPICModel::BfieldEvolveRK (
    ablastr::fields::MultiLevelVectorField const& Bfield,
    ablastr::fields::MultiLevelVectorField const& Efield,
    ablastr::fields::MultiLevelVectorField const& Jfield,
    ablastr::fields::MultiLevelScalarField const& rhofield,
    ablastr::fields::MultiLevelVectorField const& edge_lengths,
    amrex::Real dt, int lev, DtType dt_type,
    IntVect ng, std::optional<bool> nodal_sync )
{
    // Make copies of the B-field multifabs at t = n and create multifabs for
    // each direction to store the Runge-Kutta intermediate terms. Each
    // multifab has 2 components for the different terms that need to be stored.
    std::array< MultiFab, 3 > B_old;
    std::array< MultiFab, 3 > K;
    for (int ii = 0; ii < 3; ii++)
    {
        B_old[ii] = MultiFab(
            Bfield[lev][ii]->boxArray(), Bfield[lev][ii]->DistributionMap(), 1,
            Bfield[lev][ii]->nGrowVect()
        );
        MultiFab::Copy(B_old[ii], *Bfield[lev][ii], 0, 0, 1, ng);

        K[ii] = MultiFab(
            Bfield[lev][ii]->boxArray(), Bfield[lev][ii]->DistributionMap(), 2,
            Bfield[lev][ii]->nGrowVect()
        );
        K[ii].setVal(0.0);
    }

    // The Runge-Kutta scheme begins here.
    // Step 1:
    FieldPush(
        Bfield, Efield, Jfield, rhofield, edge_lengths,
        0.5_rt*dt, dt_type, ng, nodal_sync
    );

    // The Bfield is now given by:
    // B_new = B_old + 0.5 * dt * [-curl x E(B_old)] = B_old + 0.5 * dt * K0.
    for (int ii = 0; ii < 3; ii++)
    {
        // Extract 0.5 * dt * K0 for each direction into index 0 of K.
        MultiFab::LinComb(
            K[ii], 1._rt, *Bfield[lev][ii], 0, -1._rt, B_old[ii], 0, 0, 1, ng
        );
    }

    // Step 2:
    FieldPush(
        Bfield, Efield, Jfield, rhofield, edge_lengths,
        0.5_rt*dt, dt_type, ng, nodal_sync
    );

    // The Bfield is now given by:
    // B_new = B_old + 0.5 * dt * K0 + 0.5 * dt * [-curl x E(B_old + 0.5 * dt * K1)]
    //       = B_old + 0.5 * dt * K0 + 0.5 * dt * K1
    for (int ii = 0; ii < 3; ii++)
    {
        // Subtract 0.5 * dt * K0 from the Bfield for each direction, to get
        // B_new = B_old + 0.5 * dt * K1.
        MultiFab::Subtract(*Bfield[lev][ii], K[ii], 0, 0, 1, ng);
        // Extract 0.5 * dt * K1 for each direction into index 1 of K.
        MultiFab::LinComb(
            K[ii], 1._rt, *Bfield[lev][ii], 0, -1._rt, B_old[ii], 0, 1, 1, ng
        );
    }

    // Step 3:
    FieldPush(
        Bfield, Efield, Jfield, rhofield, edge_lengths,
        dt, dt_type, ng, nodal_sync
    );

    // The Bfield is now given by:
    // B_new = B_old + 0.5 * dt * K1 + dt * [-curl  x E(B_old + 0.5 * dt * K1)]
    //       = B_old + 0.5 * dt * K1 + dt * K2
    for (int ii = 0; ii < 3; ii++)
    {
        // Subtract 0.5 * dt * K1 from the Bfield for each direction to get
        // B_new = B_old + dt * K2.
        MultiFab::Subtract(*Bfield[lev][ii], K[ii], 1, 0, 1, ng);
    }

    // Step 4:
    FieldPush(
        Bfield, Efield, Jfield, rhofield, edge_lengths,
        0.5_rt*dt, dt_type, ng, nodal_sync
    );

    // The Bfield is now given by:
    // B_new = B_old + dt * K2 + 0.5 * dt * [-curl x E(B_old + dt * K2)]
    //       = B_old + dt * K2 + 0.5 * dt * K3
    for (int ii = 0; ii < 3; ii++)
    {
        // Subtract B_old from the Bfield for each direction, to get
        // B = dt * K2 + 0.5 * dt * K3.
        MultiFab::Subtract(*Bfield[lev][ii], B_old[ii], 0, 0, 1, ng);

        // Add dt * K2 + 0.5 * dt * K3 to index 0 of K (= 0.5 * dt * K0).
        MultiFab::Add(K[ii], *Bfield[lev][ii], 0, 0, 1, ng);

        // Add 2 * 0.5 * dt * K1 to index 0 of K.
        MultiFab::LinComb(
            K[ii], 1.0, K[ii], 0, 2.0, K[ii], 1, 0, 1, ng
        );

        // Overwrite the Bfield with the Runge-Kutta sum:
        // B_new = B_old + 1/3 * dt * (0.5 * K0 + K1 + K2 + 0.5 * K3).
        MultiFab::LinComb(
            *Bfield[lev][ii], 1.0, B_old[ii], 0, 1.0/3.0, K[ii], 0, 0, 1, ng
        );
    }
}

void HybridPICModel::FieldPush (
    ablastr::fields::MultiLevelVectorField const& Bfield,
    ablastr::fields::MultiLevelVectorField const& Efield,
    ablastr::fields::MultiLevelVectorField const& Jfield,
    ablastr::fields::MultiLevelScalarField const& rhofield,
    ablastr::fields::MultiLevelVectorField const& edge_lengths,
    amrex::Real dt, DtType dt_type,
    IntVect ng, std::optional<bool> nodal_sync )
{
    auto& warpx = WarpX::GetInstance();

    // Calculate J = curl x B / mu0 - J_ext
    CalculatePlasmaCurrent(Bfield, edge_lengths);
    // Calculate the E-field from Ohm's law
    HybridPICSolveE(Efield, Jfield, Bfield, rhofield, edge_lengths, true);
    warpx.FillBoundaryE(ng, nodal_sync);
    // Push forward the B-field using Faraday's law
    warpx.EvolveB(dt, dt_type);
    warpx.FillBoundaryB(ng, nodal_sync);
}
