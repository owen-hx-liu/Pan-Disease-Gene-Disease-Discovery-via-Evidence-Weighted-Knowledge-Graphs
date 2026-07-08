"""Keep Windows awake for the duration of a long compute job.

The benchmark kept dying silently mid-run because the machine went to sleep
(kernel-power event 42), which tears down the CUDA context. Calling
SetThreadExecutionState with ES_CONTINUOUS | ES_SYSTEM_REQUIRED tells Windows to
stay awake while the process lives; the state is automatically released when the
process exits, so it does NOT permanently change the user's power settings.

Usage:
    from keep_awake import keep_awake
    keep_awake()                      # stay awake until process exit
"""
import sys

ES_CONTINUOUS = 0x80000000
ES_SYSTEM_REQUIRED = 0x00000001
ES_AWAYMODE_REQUIRED = 0x00000040


def keep_awake(allow_away=False):
    """Prevent system sleep for the lifetime of this process (Windows only)."""
    if not sys.platform.startswith("win"):
        return False
    try:
        import ctypes
        flags = ES_CONTINUOUS | ES_SYSTEM_REQUIRED
        if allow_away:
            flags |= ES_AWAYMODE_REQUIRED
        ctypes.windll.kernel32.SetThreadExecutionState(flags)
        print("[keep_awake] system sleep suppressed for this run", flush=True)
        return True
    except Exception as e:
        print(f"[keep_awake] could not suppress sleep: {e}", flush=True)
        return False


def release():
    """Restore normal sleep behavior."""
    if not sys.platform.startswith("win"):
        return
    try:
        import ctypes
        ctypes.windll.kernel32.SetThreadExecutionState(ES_CONTINUOUS)
    except Exception:
        pass
