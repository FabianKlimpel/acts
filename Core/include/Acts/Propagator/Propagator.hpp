// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/algorithm/string.hpp>
#include <boost/tti/has_member_data.hpp>
#include <boost/tti/has_member_function.hpp>
#include <boost/tti/has_template.hpp>
#include <boost/tti/has_type.hpp>
#include <cmath>
#include <functional>
#include <memory>
#include <type_traits>
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/detail/LoopProtection.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Propagator/detail/VoidPropagatorComponents.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/GeometryContext.hpp"
#include "Acts/Utilities/MagneticFieldContext.hpp"
#include "Acts/Utilities/Units.hpp"

#include "Acts/Propagator/PropagatorError.hpp"
#include "Acts/Utilities/Result.hpp"

BOOST_TTI_HAS_TYPE(state_type)
BOOST_TTI_HAS_TEMPLATE(return_parameter_type, class, class)
BOOST_TTI_HAS_TYPE(State)
BOOST_TTI_HAS_MEMBER_DATA(covTransport)
BOOST_TTI_HAS_MEMBER_DATA(cov)
BOOST_TTI_HAS_MEMBER_DATA(navDir)
BOOST_TTI_HAS_MEMBER_DATA(pathAccumulated)
BOOST_TTI_HAS_MEMBER_DATA(stepSize)
BOOST_TTI_HAS_MEMBER_FUNCTION(getField)
BOOST_TTI_HAS_MEMBER_FUNCTION(position)
BOOST_TTI_HAS_MEMBER_FUNCTION(direction)
BOOST_TTI_HAS_MEMBER_FUNCTION(momentum)
BOOST_TTI_HAS_MEMBER_FUNCTION(charge)
BOOST_TTI_HAS_MEMBER_FUNCTION(surfaceReached)
BOOST_TTI_HAS_MEMBER_FUNCTION(boundState)
BOOST_TTI_HAS_MEMBER_FUNCTION(curvilinearState)
BOOST_TTI_HAS_MEMBER_FUNCTION(update)
BOOST_TTI_HAS_MEMBER_FUNCTION(covarianceTransport)
BOOST_TTI_HAS_MEMBER_FUNCTION(step)

namespace Acts {

/// @brief Simple class holding result of propagation call
///
/// @tparam parameters_t Type of final track parameters
/// @tparam result_list  Result pack for additional propagation
///                      quantities
template <typename parameters_t, typename... result_list>
struct PropagatorResult : private detail::Extendable<result_list...>
{

  /// Accessor to additional propagation quantities
  using detail::Extendable<result_list...>::get;

  /// Final track parameters - initialized to null pointer
  std::unique_ptr<const parameters_t> endParameters = nullptr;

  /// Full transport jacobian
  std::unique_ptr<const ActsMatrixD<5, 5>> transportJacobian = nullptr;

  /// Number of propagation steps that were carried out
  unsigned int steps = 0;

  /// Signed distance over which the parameters were propagated
  double pathLength = 0.;
};

/// @brief Options for propagate() call
///
/// @tparam action_list_t List of action types called after each
///    propagation step with the current propagation and stepper state
///
/// @tparam aborter_list_t List of abort conditions tested after each
///    propagation step using the current propagation and stepper state
///
template <typename action_list_t  = ActionList<>,
          typename aborter_list_t = AbortList<>>
struct PropagatorOptions
{

  /// Delete default contructor
  PropagatorOptions() = delete;

  /// PropagatorOptions copy constructor
  PropagatorOptions(const PropagatorOptions<action_list_t, aborter_list_t>& po)
      = default;

  /// PropagatorOptions with context
  PropagatorOptions(std::reference_wrapper<const GeometryContext>      gctx,
                    std::reference_wrapper<const MagneticFieldContext> mctx)
    : geoContext(gctx), magFieldContext(mctx)
  {
  }

  /// @brief Expand the Options with extended aborters
  ///
  /// @tparam extended_aborter_list_t Type of the new aborter list
  ///
  /// @param aborters The new aborter list to be used (internally)
    ///
  /// @return Propagator options with extend aborter list
  template <typename extended_aborter_list_t>
  PropagatorOptions<action_list_t, extended_aborter_list_t>
  extendAborters(extended_aborter_list_t aborters) const
  {
    PropagatorOptions<action_list_t, extended_aborter_list_t> eoptions(
        geoContext, magFieldContext);
    // Copy the options over
    eoptions.direction       = direction;
    eoptions.absPdgCode      = absPdgCode;
    eoptions.mass            = mass;
    eoptions.maxSteps        = maxSteps;
    eoptions.maxStepSize     = maxStepSize;
    eoptions.targetTolerance = targetTolerance;
    eoptions.pathLimit       = pathLimit;
    eoptions.loopProtection  = loopProtection;
    eoptions.loopFraction    = loopFraction;
    // Output option
    eoptions.debug         = debug;
    eoptions.debugString   = debugString;
    eoptions.debugPfxWidth = debugPfxWidth;
    eoptions.debugMsgWidth = debugMsgWidth;
    // Action / abort list
    eoptions.actionList = std::move(actionList);
    eoptions.abortList  = std::move(aborters);
    // And return the options
    return eoptions;
  }
  
  /// @brief Expand the Options with extended actors
  ///
  /// @tparam extended_actor_list_t Type of the new actor list
  ///
  /// @param actors The new actor list to be used (internally)
  ///
  /// @return Propagator options with extend actor list
  template <typename extended_action_list_t>
  PropagatorOptions<extended_action_list_t, aborter_list_t>
  extendActors(extended_action_list_t actors) const
  {
    PropagatorOptions<extended_action_list_t, aborter_list_t> eoptions(
        geoContext, magFieldContext);
    // Copy the options over
    eoptions.direction       = direction;
    eoptions.absPdgCode      = absPdgCode;
    eoptions.mass            = mass;
    eoptions.maxSteps        = maxSteps;
    eoptions.maxStepSize     = maxStepSize;
    eoptions.targetTolerance = targetTolerance;
    eoptions.pathLimit       = pathLimit;
    eoptions.loopProtection  = loopProtection;
    eoptions.loopFraction    = loopFraction;
    // Output option
    eoptions.debug         = debug;
    eoptions.debugString   = debugString;
    eoptions.debugPfxWidth = debugPfxWidth;
    eoptions.debugMsgWidth = debugMsgWidth;
    // Action / abort list
    eoptions.actionList = std::move(actors);
    eoptions.abortList  = std::move(abortList);
    // And return the options
    return eoptions;
  }

  /// Propagation direction
  NavigationDirection direction = forward;

  /// The |pdg| code for (eventual) material integration - pion default
  int absPdgCode = 211;

  /// The mass for the particle for (eventual) material integration
  double mass = 139.57018 * units::_MeV;

  /// Maximum number of steps for one propagate call
  unsigned int maxSteps = 1000;

  /// Absolute maximum step size
  double maxStepSize = std::numeric_limits<double>::max();

  /// Absolute maximum path length
  double pathLimit = std::numeric_limits<double>::max();

  /// Required tolerance to reach target (surface, pathlength)
  double targetTolerance = s_onSurfaceTolerance;

  /// Loop protection step, it adapts the pathLimit
  bool   loopProtection = true;
  double loopFraction   = 0.5;  ///< Allowed loop fraction, 1 is a full loop

  /// Debug output steering:
  //  -> @todo: move to a debug struct
  // - the string where debug messages are stored (optionally)
  // - it also has some formatting options
  bool        debug         = false;  ///< switch debug on
  std::string debugString   = "";     ///< the string to collect msgs
  size_t      debugPfxWidth = 30;     ///< the prefix width
  size_t      debugMsgWidth = 50;     ///< the mesage width

  // Configurations for Stepper
  /// Tolerance for the error of the integration
  double tolerance = 1e-4;

  /// Cut-off value for the step size
  double stepSizeCutOff = 0.;

  /// List of actions
  action_list_t actionList;

  /// List of abort conditions
  aborter_list_t abortList;

  /// The context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;

  /// The context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
};

/// @brief Propagator for particles (optionally in a magnetic field)
///
/// The Propagator works with a state objects given at function call
/// This state object contains the thread local state objects
///  - Navigator::state_type for object navigation and screen output
///  - Stepper::state_type state for the actual transport caching
///  (pos,dir,field)
///
/// @tparam stepper_t Type of stepper implementation of the propagation
/// @tparam naviagor_t Type of the navigator (optional)
///
/// This Propagator class serves as high-level steering code for propagating
/// track parameters. The actual implementation of the propagation has to be
/// implemented in the stepper_t object, which has to provide the following:
///
/// - a function for performing a single propagation step
/// - a type mapping for: initial track parameter type -> type of final track
///   parameters
/// - a type mapping for: (initial track parameter type and destination
///   surface type) -> type of final track parameters
/// - a type mapping for: initial track parameter type -> type of internal
///   state object
/// - a type mapping for: (initial track parameter type and destination
///   surface type) -> type of internal state object
///
template <typename stepper_t, typename navigator_t = detail::VoidNavigator>
class Propagator final
{
  using Jacobian         = ActsMatrixD<5, 5>;
  using BoundState       = std::tuple<BoundParameters, Jacobian, double>;
  using CurvilinearState = std::tuple<CurvilinearParameters, Jacobian, double>;

  // Namings
  static_assert(has_type_state_type<stepper_t>::value,
                "Stepper does not declare state_type");
  static_assert(has_template_return_parameter_type<stepper_t>::value,
                "Stepper does not declate return_parameter_type");
  // State struct
  static_assert(has_type_State<stepper_t>::value,
                "Stepper does not have a state");
  static_assert(
      has_member_data_covTransport<typename stepper_t::state_type, bool>::value,
      "StepperState does not have a covTransport variable");
  static_assert(has_member_data_cov<typename stepper_t::state_type,
                                    ActsSymMatrixD<5>>::value,
                "StepperState does not have a cov variable");
  static_assert(has_member_data_navDir<typename stepper_t::state_type,
                                       NavigationDirection>::value,
                "StepperState does not have a navDir variable");
  static_assert(has_member_data_pathAccumulated<typename stepper_t::state_type,
                                                double>::value,
                "StepperState does not have a pathAccumulated variable");
  static_assert(has_member_data_stepSize<typename stepper_t::state_type,
                                         detail::ConstrainedStep>::value,
                "StepperState does not have a stepSize variable");

  // Functions
  static_assert(
      has_member_function_getField<const stepper_t,
                                   Vector3D,
                                   boost::mpl::vector<
                                       typename stepper_t::state_type&,
                                       const Vector3D&>>::value,
      "Stepper has no getField method");
  static_assert(
      has_member_function_position<const stepper_t,
                                   Vector3D,
                                   boost::mpl::vector<const typename stepper_t::
                                                          state_type&>>::value,
      "Stepper has no position method");
  static_assert(
      has_member_function_direction<const stepper_t,
                                    Vector3D,
                                    boost::mpl::
                                        vector<const typename stepper_t::
                                                   state_type&>>::value,
      "Stepper has no direction method");
  static_assert(
      has_member_function_momentum<const stepper_t,
                                   double,
                                   boost::mpl::vector<const typename stepper_t::
                                                          state_type&>>::value,
      "Stepper has no momentum method");
  static_assert(
      has_member_function_charge<const stepper_t,
                                 double,
                                 boost::mpl::vector<const typename stepper_t::
                                                        state_type&>>::value,
      "Stepper has no charge method");
  static_assert(
      has_member_function_surfaceReached<const stepper_t,
                                         bool,
                                         boost::mpl::
                                             vector<const typename stepper_t::
                                                        state_type&,
                                                    const Surface*>>::value,
      "Stepper has no surfaceReached method");
  static_assert(
      has_member_function_boundState<const stepper_t,
                                     BoundState,
                                     boost::mpl::vector<
                                         typename stepper_t::state_type&,
                                         const Surface&,
                                         bool>>::value,
      "Stepper has no boundState method");
  static_assert(
      has_member_function_curvilinearState<const stepper_t,
                                           CurvilinearState,
                                           boost::mpl::vector<
                                               typename stepper_t::state_type&,
                                               bool>>::value,
      "Stepper has no curvilinearState method");
  static_assert(has_member_function_update<const stepper_t,
                                           void,
                                           boost::mpl::vector<
                                               typename stepper_t::state_type&,
                                               const BoundParameters&>>::value,
                "Stepper has no update method");
  static_assert(has_member_function_update<const stepper_t,
                                           void,
                                           boost::mpl::vector<
                                               typename stepper_t::state_type&,
                                               const Vector3D&,
                                               const Vector3D&,
                                               double>>::value,
                "Stepper has no update method");
  static_assert(has_member_function_covarianceTransport<const stepper_t,
                                                        void,
                                                        boost::mpl::vector<
                                                            typename stepper_t::
                                                                state_type&,
                                                            bool>>::value,
                "Stepper has no covarianceTransport method");
  static_assert(has_member_function_covarianceTransport<const stepper_t,
                                                        void,
                                                        boost::mpl::vector<
                                                            typename stepper_t::
                                                                state_type&,
                                                            const Surface&,
                                                            bool>>::value,
                "Stepper has no covarianceTransport method");

public:
  /// Type of the stepper in use for public scope
  using Stepper = stepper_t;

  /// Type of state object used by the propagation implementation
  using StepperState = typename Stepper::state_type;

  /// Typedef the navigator state
  using NavigatorState = typename navigator_t::state_type;

  /// Constructor from implementation object
  ///
  /// @param stepper The stepper implementation is moved to a private member
  /// @param navigator The navigator implementation, moved to a private member
  explicit Propagator(stepper_t stepper, navigator_t navigator = navigator_t())
    : m_stepper(std::move(stepper)), m_navigator(std::move(navigator))
  {
  }

  /// @brief private Propagator state for navigation and debugging
  ///
  /// @tparam parameters_t Type of the track parameters
  /// @tparam propagator_options_t Type of the Objections object
  ///
  /// This struct holds the common state information for propagating
  /// which is independent of the actual stepper implementation.
  template <typename propagator_options_t>
  struct State
  {

    /// Create the propagator state from the options
    ///
    /// @tparam parameters_t the type of the start parameters
    /// @tparam propagator_options_t the type of the propagator options
    ///
    /// @param start The start parameters, used to initialize stepping state
    /// @param topts The options handed over by the propagate call
    template <typename parameters_t>
    State(const parameters_t& start, const propagator_options_t& topts)
      : options(topts)
      , stepping(topts.geoContext,
                 topts.magFieldContext,
                 start,
                 topts.direction,
                 topts.maxStepSize)
      , geoContext(topts.geoContext)
    {
      // Setting the start surface
      navigation.startSurface = &start.referenceSurface();
    }

    /// These are the options - provided for each propagation step
    propagator_options_t options;

    /// Stepper state - internal state of the Stepper
    StepperState stepping;

    /// Navigation state - internal state of the Navigator
    NavigatorState navigation;

    /// Context object for the geometry
    std::reference_wrapper<const GeometryContext> geoContext;
  };

private:
  /// @brief Helper struct determining the result's type
  ///
  /// @tparam parameters_t Type of final track parameters
  /// @tparam action_list_t    List of propagation action types
  ///
  /// This helper struct provides type definitions to extract the correct
  /// propagation result type from a given TrackParameter type and an
  /// ActionList.
  ///
  template <typename parameters_t, typename action_list_t>
  struct result_type_helper
  {
    /// @brief Propagation result type for an arbitrary list of additional
    ///        propagation results
    ///
    /// @tparam args Parameter pack specifying additional propagation results
    ///
    template <typename... args>
    using this_result_type = PropagatorResult<parameters_t, args...>;

    /// @brief Propagation result type derived from a given action list
    using type = typename action_list_t::template result_type<this_result_type>;
  };

  /// @brief Short-hand type definition for propagation result derived from
  ///        an action list
  ///
  /// @tparam parameters_t Type of the final track parameters
  /// @tparam action_list_t List of propagation action types
  ///
  template <typename parameters_t, typename action_list_t>
  using action_list_t_result_t =
      typename result_type_helper<parameters_t, action_list_t>::type;

  /// @brief Propagate track parameters
  /// Private method with propagator and stepper state
  ///
  /// This function performs the propagation of the track parameters according
  /// to the internal implementation object until at least one abort condition
  /// is fulfilled, the destination surface is hit or the maximum number of
  /// steps/path length as given in the propagation options is reached.
  ///
  /// @note Does not (yet) convert into  the return_type of the propagation
  ///
  /// @tparam result_t Type of the result object for this propagation
  /// @tparam propagator_state_t Type of of propagator state with options
  ///
  /// @param [in,out] result of the propagation
  /// @param [in,out] state the propagator state object
  ///
  /// @return Propagation PropagatorStatus
  template <typename result_t, typename propagator_state_t>
  Result<result_t>
  propagate_impl(propagator_state_t& state) const
  {
    result_t result;

    // Pre-stepping call to the navigator and action list
    debugLog(state, [&] { return std::string("Entering propagation."); });

    // Navigator initialize state call
    m_navigator.status(state, m_stepper);
    // Pre-Stepping call to the action list
    state.options.actionList(state, m_stepper, result);
    // assume negative outcome, only set to true later if we actually have
    // a positive outcome.
    // This is needed for correct error logging
    bool terminatedNormally = false;
    // Pre-Stepping: abort condition check
    if (!state.options.abortList(result, state, m_stepper)) {
      // Pre-Stepping: target setting
      m_navigator.target(state, m_stepper);
      // Stepping loop
      debugLog(state, [&] { return std::string("Starting stepping loop."); });
      // Propagation loop : stepping
      for (; result.steps < state.options.maxSteps; ++result.steps) {
        // Perform a propagation step - it takes the propagation state
        Result<double> res = m_stepper.step(state);
        if (res.ok()) {
          // Accumulate the path length
          double s = *res;
          result.pathLength += s;
          debugLog(state, [&] {
            std::stringstream dstream;
            dstream << "Step with size = ";
            dstream << s;
            dstream << " performed.";
            return dstream.str();
          });
        } else {
          debugLog(state, [&] {
            std::stringstream dstream;
            dstream << "Step failed: ";
            dstream << res.error();
            return dstream.str();
          });
          // pass error to caller
          return res.error();
        }
        // Post-step
        // navigator status call - action list - aborter list - target call
        m_navigator.status(state, m_stepper);
        state.options.actionList(state, m_stepper, result);
        if (state.options.abortList(result, state, m_stepper)) {
          terminatedNormally = true;
          break;
        }
        m_navigator.target(state, m_stepper);
      }
    }

    // if we didn't terminate normally (via aborters) set navigation break.
    // this will trigger error output in the lines below
    if (!terminatedNormally) {
      state.navigation.navigationBreak = true;
    }

    // Post-stepping call to the action list
    debugLog(state, [&] { return std::string("Stepping loop done."); });
    state.options.actionList(state, m_stepper, result);

    // return progress flag here, decide on SUCCESS later
    return std::move(result);
  }

public:
  /// @brief Propagate track parameters
  ///
  /// This function performs the propagation of the track parameters using the
  /// internal stepper implementation, until at least one abort condition is
  /// fulfilled or the maximum number of steps/path length provided in the
  /// propagation options is reached.
  ///
  /// @tparam parameters_t Type of initial track parameters to propagate
  /// @tparam action_list_t Type list of actions, type ActionList<>
  /// @tparam aborter_list_t Type list of abort conditions, type AbortList<>
  /// @tparam propagator_options_t Type of the propagator options
  ///
  /// @param [in] start initial track parameters to propagate
  /// @param [in] options Propagation options, type Options<,>
  ///
  /// @return Propagation result containing the propagation status, final
  ///         track parameters, and output of actions (if they produce any)
  ///
  template <typename parameters_t,
            typename action_list_t,
            typename aborter_list_t,
            template <typename, typename> class propagator_options_t,
            typename path_aborter_t = detail::PathLimitReached>
  Result<action_list_t_result_t<
      typename stepper_t::template return_parameter_type<parameters_t>,
      action_list_t>>
  propagate(
      const parameters_t& start,
      const propagator_options_t<action_list_t, aborter_list_t>& options) const
  {

    // Type of track parameters produced by the propagation
    using ReturnParameterType =
        typename stepper_t::template return_parameter_type<parameters_t>;

    // Type of the full propagation result, including output from actions
    using ResultType
        = action_list_t_result_t<ReturnParameterType, action_list_t>;

    static_assert(std::is_copy_constructible<ReturnParameterType>::value,
                  "return track parameter type must be copy-constructible");

    // Initialize the propagation result object

    // Expand the abort list with a path aborter
    path_aborter_t pathAborter;
    auto           abortList = options.abortList.append(pathAborter);

    // The expanded options (including path limit)
    auto eOptions     = options.extendAborters(abortList);
    using OptionsType = decltype(eOptions);
    // Initialize the internal propagator state
    using StateType = State<OptionsType>;
    StateType state(start, eOptions);

    static_assert(
        has_member_function_step<const stepper_t,
                                 Result<double>,
                                 boost::mpl::vector<StateType&>>::value,
        "Step method of the Stepper is not compatible with the propagator "
        "state");

    // Apply the loop protection - it resets the internal path limit
    if (options.loopProtection) {
      detail::LoopProtection<path_aborter_t> lProtection;
      lProtection(state, m_stepper);
    }

    // Perform the actual propagation & check its outcome
    auto result = propagate_impl<ResultType>(state);

    if (result.ok()) {
      auto& propRes = *result;
      /// Convert into return type and fill the result object
      auto  curvState      = m_stepper.curvilinearState(state.stepping, true);
      auto& curvParameters = std::get<CurvilinearParameters>(curvState);
      // Fill the end parameters
      propRes.endParameters = std::make_unique<const CurvilinearParameters>(
          std::move(curvParameters));
      // Only fill the transport jacobian when covariance transport was done
      if (state.stepping.covTransport) {
        auto& tJacobian = std::get<Jacobian>(curvState);
        propRes.transportJacobian
            = std::make_unique<const Jacobian>(std::move(tJacobian));
      }
      return result;
    } else {
      return result.error();
    }
  }

  /// @brief Propagate track parameters - User method
  ///
  /// This function performs the propagation of the track parameters according
  /// to the internal implementation object until at least one abort condition
  /// is fulfilled, the destination surface is hit or the maximum number of
  /// steps/path length as given in the propagation options is reached.
  ///
  /// @tparam parameters_t Type of initial track parameters to propagate
  /// @tparam surface_t Type of target surface
  /// @tparam action_list_t Type list of actions
  /// @tparam aborter_list_t Type list of abort conditions
  /// @tparam propagator_options_t Type of the propagator options
  ///
  /// @param [in] start Initial track parameters to propagate
  /// @param [in] target Target surface of to propagate to
  /// @param [in] options Propagation options
  ///
  /// @return Propagation result containing the propagation status, final
  ///         track parameters, and output of actions (if they produce any)
  template <typename parameters_t,
            typename surface_t,
            typename action_list_t,
            typename aborter_list_t,
            template <typename, typename> class propagator_options_t,
            typename target_aborter_t = detail::SurfaceReached,
            typename path_aborter_t   = detail::PathLimitReached>
  Result<action_list_t_result_t<
      typename stepper_t::template return_parameter_type<parameters_t,
                                                         surface_t>,
      action_list_t>>
  propagate(
      const parameters_t& start,
      const surface_t&    target,
      const propagator_options_t<action_list_t, aborter_list_t>& options) const
  {

    // Type of track parameters produced at the end of the propagation
    using return_parameter_type =
        typename stepper_t::template return_parameter_type<parameters_t,
                                                           surface_t>;

    // Type of provided options
    target_aborter_t targetAborter;
    path_aborter_t   pathAborter;
    auto abortList = options.abortList.append(targetAborter, pathAborter);

    // Create the extended options and declare their type
    auto eOptions     = options.extendAborters(abortList);
    using OptionsType = decltype(eOptions);

    // Type of the full propagation result, including output from actions
    using ResultType
        = action_list_t_result_t<return_parameter_type, action_list_t>;

    // Initialize the internal propagator state
    using StateType = State<OptionsType>;
    StateType state(start, eOptions);
    state.navigation.targetSurface = &target;

    static_assert(
        has_member_function_step<const stepper_t,
                                 Result<double>,
                                 boost::mpl::vector<StateType&>>::value,
        "Step method of the Stepper is not compatible with the propagator "
        "state");

    // Apply the loop protection, it resets the interal path limit
    detail::LoopProtection<path_aborter_t> lProtection;
    lProtection(state, m_stepper);

    // Perform the actual propagation
    auto result = propagate_impl<ResultType>(state);

    if (result.ok()) {
      auto& propRes = *result;
      // Compute the final results and mark the propagation as successful
      auto  bs = m_stepper.boundState(state.stepping, target, true);
      auto& boundParameters = std::get<BoundParameters>(bs);
      // Fill the end parameters
      propRes.endParameters
          = std::make_unique<const BoundParameters>(std::move(boundParameters));
      // Only fill the transport jacobian when covariance transport was done
      if (state.stepping.covTransport) {
        auto& tJacobian = std::get<Jacobian>(bs);
        propRes.transportJacobian
            = std::make_unique<const Jacobian>(std::move(tJacobian));
      }
      return result;
    } else {
      return result.error();
    }
  }

private:
  /// Implementation of propagation algorithm
  stepper_t m_stepper;

  /// Implementation of navigator
  navigator_t m_navigator;

  /// The private propagation debug logging
  ///
  /// It needs to be fed by a lambda function that returns a string,
  /// that guarantees that the lambda is only called in the
  /// options.debug == true case in order not to spend time when not needed.
  ///
  /// @tparam propagator_state_t Type of the nested propagator state object
  ///
  /// @param state the propagator state for the debug flag, prefix/length
  /// @param logAction is a callable function that returns a stremable object
  template <typename propagator_state_t>
  void
  debugLog(propagator_state_t&                 state,
           const std::function<std::string()>& logAction) const
  {
    if (state.options.debug) {
      std::vector<std::string> lines;
      std::string              input = logAction();
      boost::split(lines, input, boost::is_any_of("\n"));
      for (const auto& line : lines) {
        std::stringstream dstream;
        dstream << "|->" << std::setw(state.options.debugPfxWidth);
        dstream << "Propagator"
                << " | ";
        dstream << std::setw(state.options.debugMsgWidth) << line << '\n';
        state.options.debugString += dstream.str();
      }
    }
  }
};

}  // namespace Acts
