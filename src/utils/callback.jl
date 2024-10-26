#* callback function for 
function get_callback_func(progress, recorder)
    default_callback_func!(state, l) = begin
        push!(recorder, (iter=state.iter, loss=l, time=now()))
        next!(progress)
        false
    end
    return default_callback_func!
end

function get_batch_callback_func(batch_size, recorder)
    cumsum_loss = 0.0
    progress = Progress(batch_size, desc="Train Epoch 1...")
    batch_callback_func!(state, l) = begin
        if state.iter % batch_size != 0
            cumsum_loss = cumsum_loss + l
            next!(progress)
        elseif state.iter % batch_size == 0
            cumsum_loss = cumsum_loss + l
            mean_loss = cumsum_loss / batch_size
            next!(progress)
            println("")
            @info Symbol(:epoch_, state.iter ÷ batch_size, Symbol(", mean_loss: "), mean_loss, ", time: ", now())
            push!(recorder, (iter=state.iter, loss=mean_loss, time=now()))
            progress = Progress(batch_size, desc="Train Epoch $(state.iter ÷ batch_size + 1)...")
            cumsum_loss = 0.0
        end
        false
    end

    return batch_callback_func!
end

function get_batch_callback_func_with_earlystop(batch_size, recorder, val_dataset, patience, min_epoch)
    cumsum_loss = 0.0
    progress = Progress(batch_size, desc="Train Epoch 1...")
    batch_callback_func!(state, l) = begin
        if state.iter % batch_size != 0
            cumsum_loss = cumsum_loss + l
            next!(progress)
        elseif state.iter % batch_size == 0
            cumsum_loss = cumsum_loss + l
            mean_loss = cumsum_loss / batch_size
            next!(progress)
            println("")
            @info Symbol(:epoch_, state.iter ÷ batch_size, Symbol(", mean_loss: "), mean_loss, ", time: ", now())
            push!(recorder, (iter=state.iter, loss=mean_loss, time=now()))
            progress = Progress(batch_size, desc="Train Epoch $(state.iter ÷ batch_size + 1)...")
            cumsum_loss = 0.0
            if state.iter ÷ batch_size >= min_epoch
                val_loss = mean([loss_func(val_dataset[key][warmup:end], tmp_pred[key][warmup:end]) for key in keys(val_dataset)])
            end
        end
        false
    end

    return batch_callback_func!
end